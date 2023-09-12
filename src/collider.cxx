// TRENTO: Reduced Thickness Event-by-event Nuclear Topology
// Copyright 2015 Jonah E. Bernhard, J. Scott Moreland
// TRENTO3D: Three-dimensional extension of TRENTO by Weiyao Ke
// MIT License

#include "collider.h"

#include <cmath>
#include <string>
#include <vector>
#include <iostream>

#include "fwd_decl.h"
#include "nucleus.h"

namespace trento {

namespace {

// Helper functions for Collider ctor.

// Create one nucleus from the configuration.
NucleusPtr create_nucleus(const VarMap& var_map, std::size_t index) {
  const auto& species = var_map["projectile"]
                        .as<std::vector<std::string>>().at(index);
  const auto& nucleon_dmin = var_map["nucleon-min-dist"].as<double>();
  return Nucleus::create(species, nucleon_dmin);
}

// Determine the maximum impact parameter.  If the configuration contains a
// non-negative value for bmax, use it; otherwise, fall back to the minimum-bias
// default.
double determine_bmax(const VarMap& var_map,
    const Nucleus& A, const Nucleus& B, const NucleonProfile& profile) {
  auto bmax = var_map["b-max"].as<double>();
  if (bmax < 0.)
    bmax = A.radius() + B.radius() + profile.max_impact();
  return bmax;
}

// Determine the asymmetry parameter (Collider::asymmetry_) for a pair of
// nuclei.  It's just rA/(rA+rB), falling back to 1/2 if both radii are zero
// (i.e. for proton-proton).
double determine_asym(const Nucleus& A, const Nucleus& B) {
  double rA = A.radius();
  double rB = B.radius();
  double sum = rA + rB;
  if (sum < 0.1)
    return 0.5;
  else
    return rA/sum;
}

}  // unnamed namespace

// Lots of members to initialize...
// Several helper functions are defined above.
Collider::Collider(const VarMap& var_map)
    : nucleon_profile_(var_map),
      output_path(var_map["output"].as<fs::path>()),
      event_(var_map),
      output_(var_map)
    //nucleusA_(create_nucleus(var_map, 0)),
    //nucleusB_(create_nucleus(var_map, 1)),
    //nevents_(var_map["number-events"].as<int>()), //dead parameter
    //bmin_(var_map["b-min"].as<double>()), //dead parameter
    //bmax_(determine_bmax(var_map, *nucleusA_, *nucleusB_, nucleon_profile_)), //dead parameter
    //asymmetry_(determine_asym(*nucleusA_, *nucleusB_)),
    {
  // Constructor body begins here.
  // Set random seed if requested.
  auto seed = var_map["random-seed"].as<int64_t>();
  if (seed > 0)
    random::engine.seed(static_cast<random::Engine::result_type>(seed));
}

// See header for explanation.
Collider::~Collider() = default;

void Collider::run_events_angantyr()
{
    fs::path nucleon_file_path = output_path.parent_path().parent_path() / "angantyr.out" / "nucleons.dat";
    if (!fs::exists(nucleon_file_path))
    {
        std::cout << "File does not exist: " << nucleon_file_path << std::endl;
        exit(1);
    }
    std::ifstream input_file(nucleon_file_path.string());

    char hash;
    int eventIndex;
    int numNucl;
    int Npart;
    double b;
    std::vector<double> target_x;
    std::vector<double> target_y;
    std::vector<double> projectile_x;
    std::vector<double> projectile_y;

    std::string line;
    getline(input_file, line);  // Read the line with number of events
    int numEvents = std::stoi(line.substr(1)); // Get the number of events

    for (int k = 0; k < numEvents; ++k) {
        getline(input_file, line);  // Read info line for each event
        std::istringstream iss(line);
        iss >> hash >> eventIndex >> numNucl >> Npart >> b;

        for (int j = 0; j < Npart; ++j) {
            getline(input_file, line);  // Read a new line for each input
            std::istringstream iss2(line);
            double x, y;
            int targetOrProjectile;
            iss2 >> x >> y >> targetOrProjectile;
            if (targetOrProjectile == 1) {
                target_x.push_back(x);
                target_y.push_back(y);
            } else {
                projectile_x.push_back(x);
                projectile_y.push_back(y);
            }
        }

        nucleusA_=Nucleus::create("#" + std::to_string(target_x.size()) ,0.0);
        nucleusB_=Nucleus::create("#" + std::to_string(projectile_x.size()) ,0.0);

        nucleusA_->sample_nucleons_angantyr(target_x, target_y);
        nucleusB_->sample_nucleons_angantyr(projectile_x, projectile_y);

        bool collision = false;
        for (auto&& A : *nucleusA_) {
            for (auto&& B : *nucleusB_) {
                collision = nucleon_profile_.participate(A, B) || collision;
            }
        }
        if(!collision)
            continue;

        event_.compute(*nucleusA_, *nucleusB_, nucleon_profile_);
        output_(k, b, event_);

        target_x.clear();
        target_y.clear();
        projectile_x.clear();
        projectile_y.clear();
    }
    input_file.close();
}


void Collider::run_events() {
  // The main event loop.
  for (int n = 0; n < nevents_; ++n) {
    // Sampling the impact parameter also implicitly prepares the nuclei for
    // event computation, i.e. by sampling nucleon positions and participants.
    double b = sample_impact_param();

    // Pass the prepared nuclei to the Event.  It computes the entropy profile
    // (thickness grid) and other event observables.
    event_.compute(*nucleusA_, *nucleusB_, nucleon_profile_);

    // Write event data.
    output_(n, b, event_);
  }
}


double Collider::sample_impact_param() {
  // Sample impact parameters until at least one nucleon-nucleon pair
  // participates.  The bool 'collision' keeps track -- it is effectively a
  // logical OR over all possible participant pairs.
  double b;
  bool collision = false;

  do {
    // Sample b from P(b)db = 2*pi*b.
    b = bmin_ + (bmax_ - bmin_) * std::sqrt(random::canonical<double>());

    // Offset each nucleus depending on the asymmetry parameter (see header).
    nucleusA_->sample_nucleons(asymmetry_ * b);
    nucleusB_->sample_nucleons((asymmetry_ - 1.) * b);

    // Check each nucleon-nucleon pair.
    for (auto&& A : *nucleusA_) {
      for (auto&& B : *nucleusB_) {
        collision = nucleon_profile_.participate(A, B) || collision;
      }
    }
  } while (!collision);

  return b;
}

}  // namespace trento
