#ifndef REGIONS_HPP
#define REGIONS_HPP

#include <memory>
#include <string>

using namespace dealii;

template <int dim>
class Region
{
public:
  Region()
  {}

  virtual ~Region() = default;

  virtual bool
  check_region(const Point<dim> &p) const = 0;
};

template <int dim>
class Tau_inclusions : public Region<dim>
{
public:
  Tau_inclusions()
  {}

  virtual bool
  check_region(const Point<dim> &p) const override
  {
    return false;
  }
};

template <>
class Tau_inclusions<3> : public Region<3>
{
public:
  Tau_inclusions()
  {}

  virtual bool
  check_region(const Point<3> &p) const override
  {
    return ((p[0] > 63.0 && p[0] < 81.0) && (p[1] > 60.0 && p[1] < 90.0) &&
            (p[2] > 46.0 && p[2] < 67.0));
  }
};

template <>
class Tau_inclusions<2> : public Region<2>
{
public:
  Tau_inclusions()
  {}

  virtual bool
  check_region(const Point<2> &p) const override
  {
    return ((p[0] > -3.0 && p[0] < 16.0) && (p[1] > -3.0 && p[1] < 10.0));
  }
};

template <int dim>
class Amyloid_Beta_deposits : public Region<dim>
{
public:
  Amyloid_Beta_deposits()
  {}

  virtual bool
  check_region(const Point<dim> &p) const override
  {
    return false;
  }
};

template <>
class Amyloid_Beta_deposits<3> : public Region<3>
{
public:
  Amyloid_Beta_deposits()
  {}

  virtual bool
  check_region(const Point<3> &p) const override
  {
    return ((p[0] > 43.0 && p[0] < 82.0) &&
            ((p[1] > 22.0 && p[1] < 90.0) || (p[1] > 100.0 && p[1] < 135.0)) &&
            (p[2] > 80.0 && p[2] < 118.0));
  }
};

template <>
class Amyloid_Beta_deposits<2> : public Region<2>
{
public:
  Amyloid_Beta_deposits()
  {}

  virtual bool
  check_region(const Point<2> &p) const override
  {
    return ((p[0] > 0.0 && p[0] < 20.0) && (p[1] > 0.0 && p[1] < 20.0));
  }
};

template <int dim>
class TPD43_inclusions : public Region<dim>
{
public:
  TPD43_inclusions()
  {}

  virtual bool
  check_region(const Point<dim> &p) const override
  {
    return false;
  }
};

template <>
class TPD43_inclusions<3> : public Region<3>
{
public:
  TPD43_inclusions()
  {}

  virtual bool
  check_region(const Point<3> &p) const override
  {
    return ((p[0] > 52.0 && p[0] < 80.0) && (p[1] > 66.0 && p[1] < 103.0) &&
            (p[2] > 80.0 && p[2] < 117.0));
  }
};

template <>
class TPD43_inclusions<2> : public Region<2>
{
public:
  TPD43_inclusions()
  {}

  virtual bool
  check_region(const Point<2> &p) const override
  {
    return ((p[0] > -15.0 && p[0] < -2.0) && (p[1] > 25.0 && p[1] < 50.0));
  }
};

template <int dim>
class Grey_matter
{
public:
  Grey_matter()
  {}

  static bool
  check_region(const Point<dim> &p)
  {
    return false;
  }
};

template <>
class Grey_matter<3>
{
public:
  Grey_matter()
  {}

  static bool
  check_region(const Point<3> &p)
  {
    return ((p[0] < 33 || p[0] > 70) || (p[1] < 25 || p[1] > 120) ||
            (p[2] > 85));
  }
};

template <>
class Grey_matter<2>
{
public:
  Grey_matter()
  {}

  static bool
  check_region(const Point<2> &p)
  {
    return (p[0] * p[0]) / (a * a) + (p[1] * p[1]) / (b * b) > 1.0;
  }

private:
  static constexpr double a = 60.0;
  static constexpr double b = 40.0;
};

template <int dim>
class Axonal_region
{
public:
  Axonal_region()
  {}

  static bool
  check_region(const Point<dim> &p)
  {
    return false;
  }

  static Tensor<1, dim>
  compute_radial_direction(const Point<dim> &p, const Point<dim> &global_center)
  {
    return Tensor<1, dim>();
  }

  static Tensor<1, dim>
  compute_circumferential_direction(const Point<dim> &p,
                                    const Point<dim> &global_center)
  {
    return Tensor<1, dim>();
  }

  static Tensor<1, dim>
  compute_axon_based_direction(const Point<dim> &p,
                               const Point<dim> &global_center)
  {
    return Tensor<1, dim>();
  }

  static Tensor<1, dim>
  get_axonal_direction(const Point<dim>  &p,
                       const Point<dim>  &global_center,
                       const std::string &orientation)
  {
    if (orientation == "radial")
      return compute_radial_direction(p, global_center);
    else if (orientation == "circumferential")
      return compute_circumferential_direction(p, global_center);
    else if (orientation == "axon-based")
      return compute_axon_based_direction(p, global_center);
    else
      throw std::invalid_argument("Invalid orientation");
  }
};

template <>
class Axonal_region<3>
{
public:
  Axonal_region()
  {}

  static bool
  check_region(const Point<3> &p)
  {
    return (p[0] < 60 && p[0] > 40) && (p[1] < 110 && p[1] > 34) &&
           (p[2] < 80 && p[2] > 50);
  }

  static Tensor<1, 3>
  compute_radial_direction(const Point<3> &p, const Point<3> &global_center)
  {
    Tensor<1, 3> radial = p - global_center;
    radial /= radial.norm();
    return radial;
  }

  static Tensor<1, 3>
  compute_circumferential_direction(const Point<3> &p,
                                    const Point<3> &global_center)
  {
    Tensor<1, 3> radial =
      Axonal_region<3>::compute_radial_direction(p, global_center);

    // Azimuthal direction perpendicular to the radial direction
    Tensor<1, 3> azimuthal;
    azimuthal[0] = -p[1];
    azimuthal[1] = p[0];
    azimuthal[2] = 0.0;
    azimuthal /= azimuthal.norm();

    // Cross product between radial and azimuthal directions
    Tensor<1, 3> circumferential;
    circumferential[0] = radial[1] * azimuthal[2] - radial[2] * azimuthal[1];
    circumferential[1] = radial[2] * azimuthal[0] - radial[0] * azimuthal[2];
    circumferential[2] = radial[0] * azimuthal[1] - radial[1] * azimuthal[0];

    return circumferential;
  }

  static Tensor<1, 3>
  compute_axon_based_direction(const Point<3> &p, const Point<3> &global_center)
  {
    if (Axonal_region<3>::check_region(p))
      {
        return compute_circumferential_direction(p, global_center);
      }
    else
      return compute_radial_direction(p, global_center);
  }

  static Tensor<1, 3>
  get_axonal_direction(const Point<3>    &p,
                       const Point<3>    &global_center,
                       const std::string &orientation)
  {
    if (orientation == "radial")
      return compute_radial_direction(p, global_center);
    else if (orientation == "circumferential")
      return compute_circumferential_direction(p, global_center);
    else if (orientation == "axon-based")
      return compute_axon_based_direction(p, global_center);
    else
      throw std::invalid_argument("Invalid orientation");
  }
};

template <>
class Axonal_region<2>
{
public:
  Axonal_region()
  {}

  static bool
  check_region(const Point<2> &p)
  {
    return (p[0] * p[0]) / (a * a) + (p[1] * p[1]) / (b * b) < 1.0;
  }

  static Tensor<1, 2>
  compute_radial_direction(const Point<2> &p, const Point<2> &global_center)
  {
    Tensor<1, 2> radial = p - global_center;
    radial /= radial.norm();
    return radial;
  }

  static Tensor<1, 2>
  compute_circumferential_direction(const Point<2> &p,
                                    const Point<2> &global_center)
  {
    Tensor<1, 2> radial =
      Axonal_region<2>::compute_radial_direction(p, global_center);

    // Azimuthal direction perpendicular to the radial direction
    Tensor<1, 2> azimuthal;
    azimuthal[0] = -p[1];
    azimuthal[1] = p[0];
    azimuthal /= azimuthal.norm();

    // Cross product between radial and azimuthal directions
    Tensor<1, 2> circumferential;
    circumferential[0] = radial[1] * azimuthal[0] - radial[0] * azimuthal[1];
    circumferential[1] = radial[0] * azimuthal[0] + radial[1] * azimuthal[1];

    return circumferential;
  }

  static Tensor<1, 2>
  compute_axon_based_direction(const Point<2> &p, const Point<2> &global_center)
  {
    if (Axonal_region<2>::check_region(p))
      {
        return compute_circumferential_direction(p, global_center);
      }
    else
      return compute_radial_direction(p, global_center);
  }

  static Tensor<1, 2>
  get_axonal_direction(const Point<2>    &p,
                       const Point<2>    &global_center,
                       const std::string &orientation)
  {
    if (orientation == "radial")
      return compute_radial_direction(p, global_center);
    else if (orientation == "circumferential")
      return compute_circumferential_direction(p, global_center);
    else if (orientation == "axon-based")
      return compute_axon_based_direction(p, global_center);
    else
      throw std::invalid_argument("Invalid orientation");
  }

private:
  static constexpr double a = 30.0;
  static constexpr double b = 20.0;
};


// Factory function to create seeding regions.
template <int dim>
std::unique_ptr<Region<dim>>
getSeedingRegion(std::string region)
{
  if (region == "Tau inclusions")
    return std::make_unique<Tau_inclusions<dim>>();
  else if (region == "Amyloid-Beta deposits")
    return std::make_unique<Amyloid_Beta_deposits<dim>>();
  else if (region == "TPD-43 inclusions")
    return std::make_unique<TPD43_inclusions<dim>>();
  else
    throw std::invalid_argument("Invalid seeding region");
}

#endif