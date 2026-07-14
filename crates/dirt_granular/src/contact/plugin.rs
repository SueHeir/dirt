use super::*;
/// Fused Hertz normal + Mindlin tangential contact force plugin.
///
/// Registers [`ContactHistoryStore`] in the [`AtomDataRegistry`] and a single
/// `hertz_mindlin_contact` system at [`ParticleSimScheduleSet::Force`].
pub struct HertzMindlinContactPlugin;

impl Plugin for HertzMindlinContactPlugin {
    fn dependencies(&self) -> Vec<std::any::TypeId> {
        grass_app::type_ids![dirt_atom::DemAtomPlugin]
    }

    fn provides(&self) -> Vec<&str> {
        vec!["contact_forces"]
    }

    fn requires(&self) -> Vec<&str> {
        vec!["dem_particles", "neighbor_list"]
    }

    fn build(&self, app: &mut App) {
        app.add_plugins(VirialStressPlugin);
        // Register ContactHistoryStore
        register_atom_data!(app, ContactHistoryStore::new());

        let contact_model = {
            let mt = app
                .get_resource_ref::<MaterialTable>()
                .expect("MaterialTable must exist before HertzMindlinContactPlugin");
            mt.contact_model.clone()
        };

        match contact_model.as_str() {
            "hooke" => {
                app.add_update_system(
                    hooke_contact_force.label(CONTACT_FORCE),
                    ParticleSimScheduleSet::Force,
                );
            }
            _ => {
                // Roadmap step 4: opt into the interior/boundary overlapped force
                // (interior pairs computed while the ghost halo is in flight). Bit-
                // identical to the standard force; set DIRT_OVERLAP_FORCE=1 to enable.
                let overlap = std::env::var("DIRT_OVERLAP_FORCE")
                    .map(|v| v != "0" && !v.is_empty())
                    .unwrap_or(false);
                if overlap {
                    app.add_update_system(
                        overlapped_contact_force.label(CONTACT_FORCE),
                        ParticleSimScheduleSet::Force,
                    );
                } else {
                    app.add_update_system(
                        hertz_mindlin_contact_force.label(CONTACT_FORCE),
                        ParticleSimScheduleSet::Force,
                    );
                }
            }
        }
    }

    fn try_build(&self, app: &mut App) -> Result<(), AppError> {
        if app.get_resource_ref::<MaterialTable>().is_none() {
            return Err(AppError::message(
                "HertzMindlinContactPlugin requires DemAtomPlugin (MaterialTable)",
            ));
        }
        self.build(app);
        Ok(())
    }
}
