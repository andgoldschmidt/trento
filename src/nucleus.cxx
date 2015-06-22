// TRENTO: Reduced Thickness Event-by-event Nuclear Topology
// Copyright 2015 Jonah E. Bernhard, J. Scott Moreland
// MIT License

#include "nucleus.h"

#include <cmath>
#include <fstream>
#include <iterator>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>

#include <boost/math/constants/constants.hpp>
#include <boost/filesystem.hpp>

#include "random.h"

namespace trento {

NucleusPtr Nucleus::create(const std::string& species) {
  // W-S params ref. in header
  // XXX: remember to add new species to the help output in main()
  if (species == "p")
    return NucleusPtr{new Proton{}};
  else if (species == "d")
    return NucleusPtr{new Deuteron{.228, 1.18}};
  else if (species == "He3")
    return NucleusPtr{new Helium3{}};
  else if (species == "Si")
    return NucleusPtr{new WoodsSaxonNucleus{28, 3.14, 0.537}};
  else if (species == "Si2")
    return NucleusPtr{new DeformedWoodsSaxonNucleus{28, 3.14, 0.537, -0.478, 0.250}};
  else if (species == "Cu")
    return NucleusPtr{new WoodsSaxonNucleus{62, 4.2, 0.596}};
  else if (species == "Cu2")
    return NucleusPtr{new DeformedWoodsSaxonNucleus{62, 4.2, 0.596, 0.162, -0.006}};
  else if (species == "Au")
    return NucleusPtr{new WoodsSaxonNucleus{197, 6.38, 0.535}};
  else if (species == "Au2")
    return NucleusPtr{new DeformedWoodsSaxonNucleus{197, 6.38, 0.535, -0.131, -0.031}};
  else if (species == "Pb")
    return NucleusPtr{new WoodsSaxonNucleus{208, 6.62, 0.546}};
  else if (species == "U")
    return NucleusPtr{new DeformedWoodsSaxonNucleus{238, 6.67, 0.44, 0.280, 0.093}};
  else
    throw std::invalid_argument{"unknown projectile species: " + species};
}

Nucleus::Nucleus(std::size_t A) : nucleons_(A), offset_(0) {}

void Nucleus::sample_nucleons(double offset) {
  offset_ = offset;
  sample_nucleons_impl();
}

void Nucleus::set_nucleon_position(Nucleon& nucleon, double x, double y) {
  nucleon.set_position(x + offset_, y);
}

Proton::Proton() : Nucleus(1) {}

/// Always zero.
double Proton::radius() const {
  return 0.;
}

/// Always place the nucleon at the origin.
void Proton::sample_nucleons_impl() {
  set_nucleon_position(*begin(), 0., 0.);
}

// ANDY{begin}
// The sample range is set at plus ten standard deviations from the
// expected value of the Hulthen distribution.  For typical values of
// a and b, the probability of sampling beyond this radius is O(10^-5).
Deuteron::Deuteron(double a, double b) : Nucleus(2), a_ (a), b_ (b),
    hulthen_dist_ (1000, 0., .25 * (1/a + 1/b + 2/(a+b))
      + 5 * std::sqrt(1/(a*a) + 1/(b*b) + 4/((a+b)*(a+b))), [a, b](double r)
      { double f = std::exp(-a*r)-std::exp(-b*r); return f*f; })
{}

/// Radius set two standard deviations past the expected value of
/// the Hulthen distribution.
double Deuteron::radius() const {
  return .5 * (1/a_ + 1/b_ + 2/(a_+b_))
    + std::sqrt(1/(a_*a_) + 1/(b_*b_) + 4/((a_+b_)*(a_+b_)));
}

/// The proton is set at (0,0) and the neutron is placed at
/// a distance determined by sampling the Hulthen distribution.
void Deuteron::sample_nucleons_impl() {
  // Place the first nucleon at the origin.
  set_nucleon_position(*begin(), 0., 0.);

  // Sample separation from Hulthen distribution.
  auto r = hulthen_dist_(random::engine);

  // Sample isotropic spherical angles.
  auto cos_theta = random::cos_theta<double>();
  auto phi = random::phi<double>();

  // Convert to transverse Cartesian coordinates
  auto r_sin_theta = r * std::sqrt(1. - cos_theta*cos_theta);
  auto x = r_sin_theta * std::cos(phi);
  auto y = r_sin_theta * std::sin(phi);

  // Place the second nucleon.
  set_nucleon_position(*std::next(begin(), 1), x, y);
  // XXX: re-center nucleon positions? seems more important now.
}

// Load file containing nucleon positions.
// Behavior is in accordance with known file structure. 
NuclearProfile::NuclearProfile(const std::string& species) {
  if (species == "He3") {
    // load file
    std::ifstream infile(he3_input_);
    if (!infile) {
      throw std::runtime_error{"file named '" + he3_input_ + "' not found."}; }
    double x1, y1, x2, y2, x3, y3, dump;
    // read in possible x,y plane positions for nucleons
    while (infile >> x1 >> y1 >> dump >> x2 >> y2 >> dump >> x3 >> y3
           >> dump >> dump >> dump >> dump) {
      positions.push_back({x1, y1, x2, y2, x3, y3}); }
  } else {
    throw std::runtime_error{"file for '" + species + "' not found."}; }
}

Helium3::Helium3() : Nucleus(3) {}

// Make nucleon profile for helium-3.
NuclearProfile Helium3::profile_ {"He3"};

// TODO: find a good default radius.
double Helium3::radius() const {
  return 0;
}

/// The three nucleons are positioned according to a random row
/// of the input file from the static object NuclearProfile.
void Helium3::sample_nucleons_impl() {
  // Grab a random row
  unsigned int row = static_cast<unsigned int>(random::canonical<>()
    * (profile_.positions.size() + 1));
  long unsigned int i = 0;
  for (auto&& nucleon : *this) {
    set_nucleon_position(nucleon, profile_.positions[row][i], 
      profile_.positions[row][i+1]);
    i+=2;
  }
}
// ANDY{end}

// Extend the W-S dist out to R + 10a; for typical values of (R, a), the
// probability of sampling a nucleon beyond this radius is O(10^-5).
WoodsSaxonNucleus::WoodsSaxonNucleus(std::size_t A, double R, double a)
    : Nucleus(A),
      R_(R),
      a_(a),
      woods_saxon_dist_(1000, 0., R + 10.*a,
        [R, a](double r) { return r*r/(1.+std::exp((r-R)/a)); })
{}

/// Return something a bit smaller than the true maximum radius.  The
/// Woods-Saxon distribution falls off very rapidly (exponentially), and since
/// this radius determines the impact parameter range, the true maximum radius
/// would cause far too many events with zero participants.
double WoodsSaxonNucleus::radius() const {
  return R_ + 3.*a_;
}

/// Sample uncorrelated Woods-Saxon nucleon positions.
void WoodsSaxonNucleus::sample_nucleons_impl() {
  for (auto&& nucleon : *this) {
    // Sample spherical radius from Woods-Saxon distribution.
    auto r = woods_saxon_dist_(random::engine);

    // Sample isotropic spherical angles.
    auto cos_theta = random::cos_theta<double>();
    auto phi = random::phi<double>();

    // Convert to transverse Cartesian coordinates
    auto r_sin_theta = r * std::sqrt(1. - cos_theta*cos_theta);
    auto x = r_sin_theta * std::cos(phi);
    auto y = r_sin_theta * std::sin(phi);

    set_nucleon_position(nucleon, x, y);
  }
  // XXX: re-center nucleon positions?
}

// Set rmax like the non-deformed case (R + 10a), but for the maximum
// "effective" radius.  The numerical coefficients for beta2 and beta4 are the
// approximate values of Y20 and Y40 at theta = 0.
DeformedWoodsSaxonNucleus::DeformedWoodsSaxonNucleus(
    std::size_t A, double R, double a, double beta2, double beta4)
    : Nucleus(A),
      R_(R),
      a_(a),
      beta2_(beta2),
      beta4_(beta4),
      rmax_(R*(1. + .63*std::fabs(beta2) + .85*std::fabs(beta4)) + 10.*a)
{}

/// Return something a bit smaller than the true maximum radius.  The
/// Woods-Saxon distribution falls off very rapidly (exponentially), and since
/// this radius determines the impact parameter range, the true maximum radius
/// would cause far too many events with zero participants.
double DeformedWoodsSaxonNucleus::radius() const {
  return rmax_ - 7.*a_;
}

double DeformedWoodsSaxonNucleus::deformed_woods_saxon_dist(
    double r, double cos_theta) const {
  auto cos_theta_sq = cos_theta*cos_theta;

  // spherical harmonics
  using math::double_constants::one_div_root_pi;
  auto Y20 = std::sqrt(5)/4. * one_div_root_pi * (3.*cos_theta_sq - 1.);
  auto Y40 = 3./16. * one_div_root_pi *
             (35.*cos_theta_sq*cos_theta_sq - 30.*cos_theta_sq + 3.);

  // "effective" radius
  auto Reff = R_ * (1. + beta2_*Y20 + beta4_*Y40);

  return 1. / (1. + std::exp((r - Reff) / a_));
}

/// Sample uncorrelated deformed Woods-Saxon nucleon positions.
void DeformedWoodsSaxonNucleus::sample_nucleons_impl() {
  // The deformed W-S distribution is defined so the symmetry axis is aligned
  // with the Z axis, so e.g. the long axis of uranium coincides with Z.
  //
  // After sampling positions, they must be randomly rotated.  In general this
  // requires three Euler rotations, but in this case we only need two
  // because there is no use in rotating about the nuclear symmetry axis.
  //
  // The two rotations are:
  //  - a polar "tilt", i.e. rotation about the X axis
  //  - an azimuthal "spin", i.e. rotation about the original Z axis

  // "tilt" angle
  const auto cos_a = random::cos_theta<double>();
  const auto sin_a = std::sqrt(1. - cos_a*cos_a);

  // "spin" angle
  const auto angle_b = random::phi<double>();
  const auto cos_b = std::cos(angle_b);
  const auto sin_b = std::sin(angle_b);

  for (auto&& nucleon : *this) {
    // Sample (r, theta) using a standard rejection method.
    // Remember to include the phase-space factors.
    double r, cos_theta;
    do {
      r = rmax_ * std::cbrt(random::canonical<double>());
      cos_theta = random::cos_theta<double>();
    } while (random::canonical<double>() > deformed_woods_saxon_dist(r, cos_theta));

    // Sample azimuthal angle.
    auto phi = random::phi<double>();

    // Convert to Cartesian coordinates.
    auto r_sin_theta = r * std::sqrt(1. - cos_theta*cos_theta);
    auto x = r_sin_theta * std::cos(phi);
    auto y = r_sin_theta * std::sin(phi);
    auto z = r * cos_theta;

    // Rotate.
    // The rotation formula was derived by composing the "tilt" and "spin"
    // rotations described above.
    auto x_rot = x*cos_b - y*cos_a*sin_b + z*sin_a*sin_b;
    auto y_rot = x*sin_b + y*cos_a*cos_b - z*sin_a*cos_b;

    set_nucleon_position(nucleon, x_rot, y_rot);
  }
}

}  // namespace trento
