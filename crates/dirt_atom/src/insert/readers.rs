use super::*;

pub(super) fn read_csv_particles(
    insert: &InsertConfig,
    file_path: &str,
    atom: &mut Atom,
    registry: &AtomDataRegistry,
    material_table: &MaterialTable,
    domain: &Domain,
    max_tag: &mut u32,
) -> Result<(), InsertFileError> {
    let mat_name = insert
        .material
        .as_deref()
        .ok_or(InsertFileError::MissingField {
            source: "CSV",
            field: "material",
        })?;
    let mat_idx = resolve_file_material(material_table, mat_name)?;

    let type_index_map = insert
        .type_map
        .as_ref()
        .map(|tm| resolve_type_map(tm, material_table))
        .transpose()?;

    // Open before checking fields that are needed only to decode rows. This
    // reports a missing input file at the fallible boundary instead of hiding
    // it behind a later configuration omission.
    let file = File::open(file_path).map_err(|e| InsertFileError::FileOpen {
        path: file_path.to_string(),
        source: e.to_string(),
    })?;

    let density = insert.density.ok_or(InsertFileError::MissingField {
        source: "CSV",
        field: "density",
    })?;

    let cols = insert.columns.clone().unwrap_or_default();
    let col_x = cols.x.unwrap_or(0);
    let col_y = cols.y.unwrap_or(1);
    let col_z = cols.z.unwrap_or(2);
    let col_radius = cols.radius;
    let col_vx = cols.vx;
    let col_vy = cols.vy;
    let col_vz = cols.vz;
    let col_atom_type = cols.atom_type;

    let default_radius = match &insert.radius {
        Some(RadiusSpec::Fixed(r)) => Some(*r),
        _ => None,
    };

    let reader = BufReader::new(file);
    let mut count = 0u32;

    for (line_num, line) in reader.lines().enumerate() {
        let line = line.map_err(|e| InsertFileError::FileRead {
            path: file_path.to_string(),
            line: line_num + 1,
            source: e.to_string(),
        })?;
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        // Skip header line if it starts with a letter
        if line_num == 0 && trimmed.chars().next().map_or(false, |c| c.is_alphabetic()) {
            continue;
        }

        let fields: Vec<&str> = trimmed.split(',').map(|s| s.trim()).collect();
        let parse = |idx: usize, name: &'static str| -> Result<f64, InsertFileError> {
            fields
                .get(idx)
                .ok_or_else(|| InsertFileError::MissingColumn {
                    path: file_path.to_string(),
                    line: line_num + 1,
                    field: name,
                })
                .and_then(|s| {
                    s.parse()
                        .map_err(|e: std::num::ParseFloatError| InsertFileError::ParseField {
                            path: file_path.to_string(),
                            line: line_num + 1,
                            field: format!("{} (column {})", name, idx),
                            value: (*s).to_string(),
                            source: e.to_string(),
                        })
                })
        };

        let x = parse(col_x, "x")?;
        let y = parse(col_y, "y")?;
        let z = parse(col_z, "z")?;
        let radius = col_radius
            .map(|c| parse(c, "radius"))
            .transpose()?
            .or(default_radius)
            .ok_or(InsertFileError::MissingDefault {
                path: file_path.to_string(),
                field: "radius",
                context: "CSV file insertion with no radius column",
            })?;
        let vx = col_vx.map(|c| parse(c, "vx")).transpose()?.unwrap_or(0.0);
        let vy = col_vy.map(|c| parse(c, "vy")).transpose()?.unwrap_or(0.0);
        let vz = col_vz.map(|c| parse(c, "vz")).transpose()?.unwrap_or(0.0);

        // Determine material: type_map lookup (if atom_type column present) → default material
        let row_mat_idx = match col_atom_type {
            Some(col) => {
                let file_type = parse(col, "atom_type")? as u32;
                lookup_material_for_type(file_type, type_index_map.as_ref(), mat_idx)
            }
            None => mat_idx,
        };
        let cutoff_padding = material_table.liquid_bridge_cutoff_padding(row_mat_idx);

        // Tag advances for every file particle (keeps tags globally consistent
        // across ranks); the atom is only stored if it lies in this subdomain.
        if owns_position(domain, &[x, y, z]) {
            insert_single_particle(
                atom,
                registry,
                DemParticle {
                    pos: [x, y, z],
                    vel: [vx, vy, vz],
                    radius,
                    cutoff_padding,
                    density,
                    mat_idx: row_mat_idx,
                    tag: *max_tag,
                },
            );
            count += 1;
        }
        *max_tag += 1;
    }

    println!(
        "DemAtomInsert: loaded {} local particles from CSV '{}'",
        count, file_path
    );
    Ok(())
}

pub(super) fn read_lammps_dump_particles(
    insert: &InsertConfig,
    file_path: &str,
    atom: &mut Atom,
    registry: &AtomDataRegistry,
    material_table: &MaterialTable,
    domain: &Domain,
    max_tag: &mut u32,
) -> Result<(), InsertFileError> {
    let mat_name = insert
        .material
        .as_deref()
        .ok_or(InsertFileError::MissingField {
            source: "lammps_dump",
            field: "material",
        })?;
    let mat_idx = resolve_file_material(material_table, mat_name)?;

    let type_index_map = insert
        .type_map
        .as_ref()
        .map(|tm| resolve_type_map(tm, material_table))
        .transpose()?;

    let density = insert.density.ok_or(InsertFileError::MissingField {
        source: "lammps_dump",
        field: "density",
    })?;

    let default_radius = match &insert.radius {
        Some(RadiusSpec::Fixed(r)) => Some(*r),
        _ => None,
    };

    let file = File::open(file_path).map_err(|e| InsertFileError::FileOpen {
        path: file_path.to_string(),
        source: e.to_string(),
    })?;
    let reader = BufReader::new(file);
    let mut lines = reader.lines();

    // Parse LAMMPS dump format
    let mut n_atoms: usize = 0;
    let mut column_names: Vec<String> = Vec::new();
    let mut reading_atoms = false;
    let mut count = 0u32;

    // Helper to find column index by name
    let find_col =
        |names: &[String], name: &str| -> Option<usize> { names.iter().position(|n| n == name) };

    let mut line_num = 0usize;
    while let Some(line) = lines.next() {
        line_num += 1;
        let line = line.map_err(|e| InsertFileError::FileRead {
            path: file_path.to_string(),
            line: line_num,
            source: e.to_string(),
        })?;
        let trimmed = line.trim();

        if trimmed == "ITEM: NUMBER OF ATOMS" {
            if let Some(next) = lines.next() {
                line_num += 1;
                let next = next.map_err(|e| InsertFileError::FileRead {
                    path: file_path.to_string(),
                    line: line_num,
                    source: e.to_string(),
                })?;
                n_atoms = next.trim().parse().map_err(|e: std::num::ParseIntError| {
                    InsertFileError::ParseField {
                        path: file_path.to_string(),
                        line: line_num,
                        field: "number of atoms".to_string(),
                        value: next.trim().to_string(),
                        source: e.to_string(),
                    }
                })?;
            }
            continue;
        }

        if trimmed.starts_with("ITEM: ATOMS") {
            // Parse column names from header: "ITEM: ATOMS id type x y z ..."
            column_names = trimmed
                .strip_prefix("ITEM: ATOMS")
                .unwrap_or("")
                .split_whitespace()
                .map(|s| s.to_string())
                .collect();
            reading_atoms = true;
            continue;
        }

        if trimmed.starts_with("ITEM:") {
            reading_atoms = false;
            continue;
        }

        if reading_atoms && !trimmed.is_empty() {
            let fields: Vec<&str> = trimmed.split_whitespace().collect();
            if fields.len() < column_names.len() {
                continue;
            }

            let parse_col = |name: &'static str| -> Result<Option<f64>, InsertFileError> {
                let Some(i) = find_col(&column_names, name) else {
                    return Ok(None);
                };
                let Some(value) = fields.get(i) else {
                    return Err(InsertFileError::MissingColumn {
                        path: file_path.to_string(),
                        line: line_num,
                        field: name,
                    });
                };
                value
                    .parse()
                    .map(Some)
                    .map_err(|e: std::num::ParseFloatError| InsertFileError::ParseField {
                        path: file_path.to_string(),
                        line: line_num,
                        field: name.to_string(),
                        value: (*value).to_string(),
                        source: e.to_string(),
                    })
            };

            let x = parse_col("x")?.ok_or(InsertFileError::MissingColumn {
                path: file_path.to_string(),
                line: line_num,
                field: "x",
            })?;
            let y = parse_col("y")?.ok_or(InsertFileError::MissingColumn {
                path: file_path.to_string(),
                line: line_num,
                field: "y",
            })?;
            let z = parse_col("z")?.ok_or(InsertFileError::MissingColumn {
                path: file_path.to_string(),
                line: line_num,
                field: "z",
            })?;
            let vx = parse_col("vx")?.unwrap_or(0.0);
            let vy = parse_col("vy")?.unwrap_or(0.0);
            let vz = parse_col("vz")?.unwrap_or(0.0);
            let radius =
                parse_col("radius")?
                    .or(default_radius)
                    .ok_or(InsertFileError::MissingDefault {
                        path: file_path.to_string(),
                        field: "radius",
                        context: "LAMMPS dump file insertion with no radius column",
                    })?;

            // Determine material: type_map override → default material
            let row_mat_idx = match parse_col("type")? {
                Some(t) => lookup_material_for_type(t as u32, type_index_map.as_ref(), mat_idx),
                None => mat_idx,
            };
            let cutoff_padding = material_table.liquid_bridge_cutoff_padding(row_mat_idx);

            if owns_position(domain, &[x, y, z]) {
                insert_single_particle(
                    atom,
                    registry,
                    DemParticle {
                        pos: [x, y, z],
                        vel: [vx, vy, vz],
                        radius,
                        cutoff_padding,
                        density,
                        mat_idx: row_mat_idx,
                        tag: *max_tag,
                    },
                );
                count += 1;
            }
            *max_tag += 1;
        }
    }

    let _ = n_atoms; // used for format validation if needed
    println!(
        "DemAtomInsert: loaded {} local particles from LAMMPS dump '{}'",
        count, file_path
    );
    Ok(())
}

/// Parse a field from a LAMMPS data file, with a user-friendly error on failure.
pub(super) fn parse_field<T: std::str::FromStr>(
    value: &str,
    field_name: &str,
    line_num: usize,
    file_path: &str,
) -> Result<T, InsertFileError>
where
    T::Err: std::fmt::Display,
{
    value.parse::<T>().map_err(|e| InsertFileError::ParseField {
        path: file_path.to_string(),
        line: line_num,
        field: field_name.to_string(),
        value: value.to_string(),
        source: e.to_string(),
    })
}

pub(super) fn read_lammps_data_particles(
    insert: &InsertConfig,
    file_path: &str,
    atom: &mut Atom,
    registry: &AtomDataRegistry,
    material_table: &MaterialTable,
    domain: &Domain,
    max_tag: &mut u32,
) -> Result<(), InsertFileError> {
    let mat_name = insert
        .material
        .as_deref()
        .ok_or(InsertFileError::MissingField {
            source: "lammps_data",
            field: "material",
        })?;
    let mat_idx = resolve_file_material(material_table, mat_name)?;

    let type_index_map = insert
        .type_map
        .as_ref()
        .map(|tm| resolve_type_map(tm, material_table))
        .transpose()?;

    let default_density = insert.density;
    let default_radius = match &insert.radius {
        Some(RadiusSpec::Fixed(r)) => Some(*r),
        _ => None,
    };

    let file = File::open(file_path).map_err(|e| InsertFileError::FileOpen {
        path: file_path.to_string(),
        source: e.to_string(),
    })?;
    let reader = BufReader::new(file);
    let lines: Vec<String> = reader
        .lines()
        .enumerate()
        .map(|(i, l)| {
            l.map_err(|e| InsertFileError::FileRead {
                path: file_path.to_string(),
                line: i + 1,
                source: e.to_string(),
            })
        })
        .collect::<Result<_, _>>()?;

    // Detect atom style from config or from "Atoms # style" header
    let config_style = insert.atom_style.as_deref();

    // Find section start indices
    let mut atoms_start = None;
    let mut atoms_style = None;
    let mut velocities_start = None;

    for (i, line) in lines.iter().enumerate() {
        let trimmed = line.trim();
        if trimmed.starts_with("Atoms") {
            atoms_start = Some(i + 1);
            // Try to detect style from "Atoms # style" comment
            if let Some(comment) = trimmed.strip_prefix("Atoms") {
                let comment = comment.trim();
                if let Some(style) = comment.strip_prefix('#') {
                    let style = style.trim();
                    if !style.is_empty() {
                        atoms_style = Some(style.to_string());
                    }
                }
            }
        } else if trimmed == "Velocities" {
            velocities_start = Some(i + 1);
        }
    }

    let atom_style = config_style
        .map(|s| s.to_string())
        .or(atoms_style)
        .unwrap_or_else(|| "atomic".to_string());

    let atoms_start = atoms_start.ok_or(InsertFileError::MissingSection {
        path: file_path.to_string(),
        section: "Atoms",
    })?;

    // Parse Atoms section
    struct ParsedAtom {
        id: u32,
        atom_type: u32,
        pos: [f64; 3],
        radius: f64,
        density: f64,
    }

    let section_headers = [
        "Atoms",
        "Velocities",
        "Bonds",
        "Angles",
        "Dihedrals",
        "Impropers",
        "Masses",
        "Pair Coeffs",
    ];
    let is_section_header = |line: &str| -> bool {
        let trimmed = line.trim();
        section_headers.iter().any(|h| trimmed.starts_with(h))
    };

    let mut parsed_atoms: Vec<ParsedAtom> = Vec::new();

    for i in atoms_start..lines.len() {
        let trimmed = lines[i].trim();
        if trimmed.is_empty() {
            continue;
        }
        if is_section_header(trimmed) {
            break;
        }
        // Skip comment lines
        if trimmed.starts_with('#') {
            continue;
        }

        let fields: Vec<&str> = trimmed.split_whitespace().collect();

        match atom_style.as_str() {
            "atomic" => {
                // id type x y z
                if fields.len() < 5 {
                    return Err(InsertFileError::RowTooShort {
                        path: file_path.to_string(),
                        line: i + 1,
                        style: "atomic".to_string(),
                        expected: 5,
                        found: fields.len(),
                    });
                }
                let id: u32 = parse_field(fields[0], "atom id", i + 1, file_path)?;
                let atype: u32 = parse_field(fields[1], "atom type", i + 1, file_path)?;
                let x: f64 = parse_field(fields[2], "x coordinate", i + 1, file_path)?;
                let y: f64 = parse_field(fields[3], "y coordinate", i + 1, file_path)?;
                let z: f64 = parse_field(fields[4], "z coordinate", i + 1, file_path)?;
                let radius = default_radius.ok_or(InsertFileError::MissingDefault {
                    path: file_path.to_string(),
                    field: "radius",
                    context: "atomic style LAMMPS data",
                })?;
                let density = default_density.ok_or(InsertFileError::MissingDefault {
                    path: file_path.to_string(),
                    field: "density",
                    context: "atomic style LAMMPS data",
                })?;
                parsed_atoms.push(ParsedAtom {
                    id,
                    atom_type: atype,
                    pos: [x, y, z],
                    radius,
                    density,
                });
            }
            "sphere" | "bpm/sphere" => {
                // id type diameter density x y z
                if fields.len() < 7 {
                    return Err(InsertFileError::RowTooShort {
                        path: file_path.to_string(),
                        line: i + 1,
                        style: atom_style.clone(),
                        expected: 7,
                        found: fields.len(),
                    });
                }
                let id: u32 = parse_field(fields[0], "atom id", i + 1, file_path)?;
                let atype: u32 = parse_field(fields[1], "atom type", i + 1, file_path)?;
                let diameter: f64 = parse_field(fields[2], "diameter", i + 1, file_path)?;
                let density: f64 = parse_field(fields[3], "density", i + 1, file_path)?;
                let x: f64 = parse_field(fields[4], "x coordinate", i + 1, file_path)?;
                let y: f64 = parse_field(fields[5], "y coordinate", i + 1, file_path)?;
                let z: f64 = parse_field(fields[6], "z coordinate", i + 1, file_path)?;
                parsed_atoms.push(ParsedAtom {
                    id,
                    atom_type: atype,
                    pos: [x, y, z],
                    radius: diameter / 2.0,
                    density,
                });
            }
            other => {
                return Err(InsertFileError::UnsupportedAtomStyle {
                    path: file_path.to_string(),
                    style: other.to_string(),
                });
            }
        }
    }

    // Parse Velocities section (optional) — build id → [vx, vy, vz] map
    let mut velocity_map: HashMap<u32, [f64; 3]> = HashMap::new();
    if let Some(vel_start) = velocities_start {
        for i in vel_start..lines.len() {
            let trimmed = lines[i].trim();
            if trimmed.is_empty() {
                continue;
            }
            if is_section_header(trimmed) {
                break;
            }
            if trimmed.starts_with('#') {
                continue;
            }
            let fields: Vec<&str> = trimmed.split_whitespace().collect();
            if fields.len() >= 4 {
                let id: u32 = parse_field(fields[0], "atom id (Velocities)", i + 1, file_path)?;
                let vx: f64 = parse_field(fields[1], "vx", i + 1, file_path)?;
                let vy: f64 = parse_field(fields[2], "vy", i + 1, file_path)?;
                let vz: f64 = parse_field(fields[3], "vz", i + 1, file_path)?;
                velocity_map.insert(id, [vx, vy, vz]);
            }
        }
    }

    // Insert all parsed atoms (only those owned by this subdomain).
    let mut count = 0usize;
    for pa in parsed_atoms {
        let vel = velocity_map.get(&pa.id).copied().unwrap_or([0.0; 3]);
        let row_mat_idx = lookup_material_for_type(pa.atom_type, type_index_map.as_ref(), mat_idx);
        let cutoff_padding = material_table.liquid_bridge_cutoff_padding(row_mat_idx);
        if owns_position(domain, &pa.pos) {
            insert_single_particle(
                atom,
                registry,
                DemParticle {
                    pos: pa.pos,
                    vel,
                    radius: pa.radius,
                    cutoff_padding,
                    density: pa.density,
                    mat_idx: row_mat_idx,
                    tag: *max_tag,
                },
            );
            count += 1;
        }
        *max_tag += 1;
    }

    println!(
        "DemAtomInsert: loaded {} local particles from LAMMPS data file '{}' (style: {})",
        count, file_path, atom_style
    );
    Ok(())
}
