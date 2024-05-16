#ifndef REGIONS_HPP
#define REGIONS_HPP

#include <memory>
#include <string>


using namespace dealii;


template <int dim>
class Region
{
public:
  Region()          = default;
  virtual ~Region() = default;

  virtual bool
  check_region(const Point<dim> &p) const = 0;
};

template <int dim>
class Tau_inclusions : public Region<dim>
{
public:
  Tau_inclusions() = default;

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
  Tau_inclusions() = default;

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
  Tau_inclusions() = default;

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
  Amyloid_Beta_deposits() = default;

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
  Amyloid_Beta_deposits() = default;

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
  Amyloid_Beta_deposits() = default;

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
  TPD43_inclusions() = default;

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
  TPD43_inclusions() = default;

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
  TPD43_inclusions() = default;

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
  Grey_matter() = default;

  static bool
  check_region(const Point<dim> &center)
  {
    return false;
  }
};

template <>
class Grey_matter<3>
{
public:
  Grey_matter() = default;

  static bool
  check_region(const Point<3> &center)
  {
    return ((center[0] < 33 || center[0] > 70) ||
            (center[1] < 25 || center[1] > 120) || (center[2] > 85));
  }
};

template <>
class Grey_matter<2>
{
public:
  Grey_matter() = default;

  static bool
  check_region(const Point<2> &center)
  {
    return ((center[0] < 33 || center[0] > 70) ||
            (center[1] < 25 || center[1] > 120));
  }
};

template <int dim>
class Axonal_region
{
public:
  Axonal_region()
  {}

  static bool
  check_region(const Point<dim> &center)
  {
    return false;
  }
};

template <>
class Axonal_region<3>
{
public:
  Axonal_region()
  {}

  static bool
  check_region(const Point<3> &center)
  {
    return (center[0] < 60 && center[0] > 40) &&
           (center[1] < 110 && center[1] > 34) &&
           (center[2] < 80 && center[2] > 50);
  }
};

template <>
class Axonal_region<2>
{
public:
  Axonal_region()
  {}

  static bool
  check_region(const Point<2> &center)
  {
    return (center[0] < 60 && center[0] > 40) &&
           (center[1] < 110 && center[1] > 34);
  }
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