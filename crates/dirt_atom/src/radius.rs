//! Radius specifications: fixed or distribution-based particle radius.

use rand::Rng;
use rand_distr::{Distribution, LogNormal, Normal};
use serde::Deserialize;

// ── RadiusSpec — fixed or distribution-based particle radius ─────────────

/// Particle radius specification: either a fixed value or a statistical distribution.
///
/// In TOML, use a plain number for fixed radius or a table with `distribution` key:
/// ```toml
/// radius = 0.001
/// radius = { distribution = "uniform", min = 0.0008, max = 0.0012 }
/// ```
#[derive(Deserialize, Clone, Debug)]
#[serde(untagged)]
pub enum RadiusSpec {
    /// A single fixed radius (m) shared by every inserted particle.
    Fixed(f64),
    /// A statistical distribution sampled per particle.
    Distribution(RadiusDistribution),
}

/// Statistical distribution for particle radii.
///
/// For `lognormal`, `mean` and `std` are the desired mean and standard deviation
/// of the actual radius distribution (not the underlying normal parameters).
#[derive(Deserialize, Clone, Debug)]
#[serde(tag = "distribution", rename_all = "lowercase")]
pub enum RadiusDistribution {
    /// Uniform distribution on `[min, max]` (m).
    Uniform {
        /// Inclusive lower bound of the radius range (m).
        min: f64,
        /// Inclusive upper bound of the radius range (m).
        max: f64,
    },
    /// Normal distribution with the given mean and standard deviation (m).
    Gaussian {
        /// Mean radius (m).
        mean: f64,
        /// Standard deviation of the radius (m).
        std: f64,
    },
    /// Log-normal distribution whose realized radii have the given mean and std (m).
    Lognormal {
        /// Desired mean of the sampled radii (m).
        mean: f64,
        /// Desired standard deviation of the sampled radii (m).
        std: f64,
    },
    /// Discrete set of radii sampled according to the given weights.
    Discrete {
        /// Candidate radii (m).
        values: Vec<f64>,
        /// Relative sampling weights, aligned with `values`.
        weights: Vec<f64>,
    },
}

impl RadiusSpec {
    /// Sample a radius from this specification, returning a config-scoped error
    /// instead of panicking when user-provided distribution parameters are invalid.
    pub fn try_sample(&self, rng: &mut impl Rng) -> Result<f64, String> {
        match self {
            RadiusSpec::Fixed(r) => Ok(*r),
            RadiusSpec::Distribution(d) => d.try_sample(rng),
        }
    }

    /// Sample a radius from this specification.
    pub fn sample(&self, rng: &mut impl Rng) -> f64 {
        // Internal/test convenience wrapper. Runtime config paths use
        // `try_sample` so malformed user distributions become typed errors.
        self.try_sample(rng)
            .expect("radius distribution parameters must be validated before sampling")
    }

    /// Conservative upper bound on radius (for spatial hash cell sizing).
    pub fn max_radius(&self) -> f64 {
        match self {
            RadiusSpec::Fixed(r) => *r,
            RadiusSpec::Distribution(d) => d.max_radius(),
        }
    }
}

impl RadiusDistribution {
    fn try_sample(&self, rng: &mut impl Rng) -> Result<f64, String> {
        match self {
            RadiusDistribution::Uniform { min, max } => {
                if min >= max {
                    return Err(format!(
                        "uniform radius requires min < max, got min={} max={}",
                        min, max
                    ));
                }
                Ok(rng.random_range(*min..*max))
            }
            RadiusDistribution::Gaussian { mean, std } => {
                let normal = Normal::new(*mean, *std).map_err(|e| {
                    format!("invalid Gaussian radius parameters (mean={mean}, std={std}): {e}")
                })?;
                Ok(normal.sample(rng).max(1e-15)) // clamp to positive
            }
            RadiusDistribution::Lognormal { mean, std } => {
                if *mean <= 0.0 || *std <= 0.0 {
                    return Err(format!(
                        "lognormal radius requires mean and std > 0, got mean={} std={}",
                        mean, std
                    ));
                }
                // Convert actual mean/std to underlying normal parameters
                let sigma_sq = (1.0 + (std / mean).powi(2)).ln();
                let mu = mean.ln() - sigma_sq / 2.0;
                let sigma = sigma_sq.sqrt();
                let ln = LogNormal::new(mu, sigma).map_err(|e| {
                    format!("invalid lognormal radius parameters (mean={mean}, std={std}): {e}")
                })?;
                Ok(ln.sample(rng))
            }
            RadiusDistribution::Discrete { values, weights } => {
                if values.is_empty() {
                    return Err("discrete radius requires at least one value".to_string());
                }
                if values.len() != weights.len() {
                    return Err(format!(
                        "discrete radius requires values/weights length match, got {} values and {} weights",
                        values.len(),
                        weights.len()
                    ));
                }
                let total: f64 = weights.iter().sum();
                if total <= 0.0 {
                    return Err(format!(
                        "discrete radius requires positive total weight, got {}",
                        total
                    ));
                }
                let r: f64 = rng.random_range(0.0..total);
                let mut cumulative = 0.0;
                for (i, w) in weights.iter().enumerate() {
                    cumulative += w;
                    if r < cumulative {
                        return Ok(values[i]);
                    }
                }
                // Non-empty was checked above; this only handles roundoff at the
                // upper edge of the cumulative distribution.
                Ok(*values
                    .last()
                    .expect("discrete distribution was validated as non-empty"))
            }
        }
    }

    fn max_radius(&self) -> f64 {
        match self {
            RadiusDistribution::Uniform { max, .. } => *max,
            RadiusDistribution::Gaussian { mean, std } => mean + 4.0 * std,
            RadiusDistribution::Lognormal { mean, std } => mean + 4.0 * std,
            RadiusDistribution::Discrete { values, .. } => {
                values.iter().cloned().fold(f64::NEG_INFINITY, f64::max)
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use soil_core::toml;

    #[test]
    fn radius_spec_fixed_deserialization() {
        let toml_str = "radius = 0.001";
        #[derive(Deserialize)]
        struct Wrapper {
            radius: RadiusSpec,
        }
        let w: Wrapper = toml::from_str(toml_str).unwrap();
        match w.radius {
            RadiusSpec::Fixed(r) => assert!((r - 0.001).abs() < 1e-15),
            _ => panic!("Expected Fixed variant"),
        }
    }

    #[test]
    fn radius_spec_uniform_deserialization() {
        let toml_str = r#"radius = { distribution = "uniform", min = 0.0008, max = 0.0012 }"#;
        #[derive(Deserialize)]
        struct Wrapper {
            radius: RadiusSpec,
        }
        let w: Wrapper = toml::from_str(toml_str).unwrap();
        match &w.radius {
            RadiusSpec::Distribution(RadiusDistribution::Uniform { min, max }) => {
                assert!((min - 0.0008).abs() < 1e-15);
                assert!((max - 0.0012).abs() < 1e-15);
            }
            other => panic!("Expected Uniform, got {:?}", other),
        }
    }

    #[test]
    fn radius_spec_gaussian_deserialization() {
        let toml_str = r#"radius = { distribution = "gaussian", mean = 0.001, std = 0.0001 }"#;
        #[derive(Deserialize)]
        struct Wrapper {
            radius: RadiusSpec,
        }
        let w: Wrapper = toml::from_str(toml_str).unwrap();
        match &w.radius {
            RadiusSpec::Distribution(RadiusDistribution::Gaussian { mean, std }) => {
                assert!((mean - 0.001).abs() < 1e-15);
                assert!((std - 0.0001).abs() < 1e-15);
            }
            other => panic!("Expected Gaussian, got {:?}", other),
        }
    }

    #[test]
    fn radius_spec_lognormal_deserialization() {
        let toml_str = r#"radius = { distribution = "lognormal", mean = 0.001, std = 0.0001 }"#;
        #[derive(Deserialize)]
        struct Wrapper {
            radius: RadiusSpec,
        }
        let w: Wrapper = toml::from_str(toml_str).unwrap();
        match &w.radius {
            RadiusSpec::Distribution(RadiusDistribution::Lognormal { mean, std }) => {
                assert!((mean - 0.001).abs() < 1e-15);
                assert!((std - 0.0001).abs() < 1e-15);
            }
            other => panic!("Expected Lognormal, got {:?}", other),
        }
    }

    #[test]
    fn radius_spec_discrete_deserialization() {
        let toml_str = r#"radius = { distribution = "discrete", values = [0.001, 0.0015], weights = [0.7, 0.3] }"#;
        #[derive(Deserialize)]
        struct Wrapper {
            radius: RadiusSpec,
        }
        let w: Wrapper = toml::from_str(toml_str).unwrap();
        match &w.radius {
            RadiusSpec::Distribution(RadiusDistribution::Discrete { values, weights }) => {
                assert_eq!(values.len(), 2);
                assert_eq!(weights.len(), 2);
                assert!((values[0] - 0.001).abs() < 1e-15);
                assert!((weights[0] - 0.7).abs() < 1e-15);
            }
            other => panic!("Expected Discrete, got {:?}", other),
        }
    }

    #[test]
    fn radius_spec_sampling_fixed() {
        let spec = RadiusSpec::Fixed(0.005);
        let mut rng = rand::rng();
        for _ in 0..10 {
            assert!((spec.sample(&mut rng) - 0.005).abs() < 1e-15);
        }
    }

    #[test]
    fn radius_spec_sampling_uniform() {
        let spec = RadiusSpec::Distribution(RadiusDistribution::Uniform {
            min: 0.001,
            max: 0.002,
        });
        let mut rng = rand::rng();
        for _ in 0..100 {
            let r = spec.sample(&mut rng);
            assert!(r >= 0.001 && r < 0.002, "uniform sample {} out of range", r);
        }
    }

    #[test]
    fn radius_spec_sampling_gaussian() {
        let spec = RadiusSpec::Distribution(RadiusDistribution::Gaussian {
            mean: 0.01,
            std: 0.001,
        });
        let mut rng = rand::rng();
        let samples: Vec<f64> = (0..1000).map(|_| spec.sample(&mut rng)).collect();
        let mean: f64 = samples.iter().sum::<f64>() / samples.len() as f64;
        assert!(
            (mean - 0.01).abs() < 0.001,
            "gaussian mean should be ~0.01, got {}",
            mean
        );
    }

    #[test]
    fn radius_spec_sampling_lognormal() {
        let spec = RadiusSpec::Distribution(RadiusDistribution::Lognormal {
            mean: 0.01,
            std: 0.001,
        });
        let mut rng = rand::rng();
        let samples: Vec<f64> = (0..5000).map(|_| spec.sample(&mut rng)).collect();
        let mean: f64 = samples.iter().sum::<f64>() / samples.len() as f64;
        // Lognormal mean should match the requested mean
        assert!(
            (mean - 0.01).abs() < 0.002,
            "lognormal mean should be ~0.01, got {}",
            mean
        );
        // All samples should be positive
        assert!(
            samples.iter().all(|&r| r > 0.0),
            "lognormal samples should all be positive"
        );
    }

    #[test]
    fn radius_spec_sampling_discrete() {
        let spec = RadiusSpec::Distribution(RadiusDistribution::Discrete {
            values: vec![0.001, 0.002],
            weights: vec![0.7, 0.3],
        });
        let mut rng = rand::rng();
        let mut count_small = 0;
        let n = 10000;
        for _ in 0..n {
            let r = spec.sample(&mut rng);
            assert!(
                (r - 0.001).abs() < 1e-15 || (r - 0.002).abs() < 1e-15,
                "discrete sample should be one of the values"
            );
            if (r - 0.001).abs() < 1e-15 {
                count_small += 1;
            }
        }
        let ratio = count_small as f64 / n as f64;
        assert!(
            (ratio - 0.7).abs() < 0.05,
            "discrete ratio should be ~0.7, got {}",
            ratio
        );
    }

    #[test]
    fn malformed_radius_distribution_reports_error() {
        let mut rng = rand::rng();
        let spec = RadiusSpec::Distribution(RadiusDistribution::Lognormal {
            mean: 0.0,
            std: 0.1,
        });

        let err = spec
            .try_sample(&mut rng)
            .expect_err("invalid lognormal config should not panic");
        assert!(err.contains("lognormal radius requires mean and std > 0"));
    }

    #[test]
    fn empty_discrete_radius_distribution_reports_error() {
        let mut rng = rand::rng();
        let spec = RadiusSpec::Distribution(RadiusDistribution::Discrete {
            values: vec![],
            weights: vec![],
        });

        let err = spec
            .try_sample(&mut rng)
            .expect_err("empty discrete config should not panic");
        assert!(err.contains("discrete radius requires at least one value"));
    }

    #[test]
    fn radius_spec_max_radius() {
        assert!((RadiusSpec::Fixed(0.005).max_radius() - 0.005).abs() < 1e-15);
        assert!(
            (RadiusSpec::Distribution(RadiusDistribution::Uniform {
                min: 0.001,
                max: 0.003
            })
            .max_radius()
                - 0.003)
                .abs()
                < 1e-15
        );
        assert!(
            (RadiusSpec::Distribution(RadiusDistribution::Discrete {
                values: vec![0.001, 0.005, 0.002],
                weights: vec![1.0, 1.0, 1.0],
            })
            .max_radius()
                - 0.005)
                .abs()
                < 1e-15
        );
    }
}
