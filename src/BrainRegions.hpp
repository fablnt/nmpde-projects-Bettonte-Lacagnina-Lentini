#ifndef REGIONS_HPP
#define REGIONS_HPP

#include <memory>
#include <string>

using namespace dealii;
/**
 * Base class for seeding regions in the brain.
 *
 * @tparam dim Dimension of the problem.
 */
template <int dim>
class SeedingRegion
{
public:
  SeedingRegion()
  {}

  virtual ~SeedingRegion() = default;

  /**
   * Check if a point is inside the seeding region.
   *
   * @param p Point to check.
   * @return True if the point is inside the region, false otherwise.
   *
   * @note Pure virtual function that must be implemented by derived classes.
   */
  virtual bool
  check_region(const Point<dim> &p) const = 0;
};

/**
 * Derived class for Tau inclusions in the brain.
 */
template <int dim>
class Tau_inclusions : public SeedingRegion<dim>
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

/**
 * Specialization of the Tau_inclusions class for 3D problems.
 */
template <>
class Tau_inclusions<3> : public SeedingRegion<3>
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

/**
 * Specialization of the Tau_inclusions class for 2D problems.
 */
template <>
class Tau_inclusions<2> : public SeedingRegion<2>
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

/**
 * Derived class for Amyloid Beta deposits in the brain.
 */
template <int dim>
class Amyloid_Beta_deposits : public SeedingRegion<dim>
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

/**
 * Specialization of the Amyloid_Beta_deposits class for 3D problems.
 */
template <>
class Amyloid_Beta_deposits<3> : public SeedingRegion<3>
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

/**
 * Specialization of the Amyloid_Beta_deposits class for 2D problems.
 */
template <>
class Amyloid_Beta_deposits<2> : public SeedingRegion<2>
{
public:
  Amyloid_Beta_deposits()
  {}

  virtual bool
  check_region(const Point<2> &p) const override
  {
    return ((p[0] > -70 && p[0] < -12) ||
            ((p[0] > -12 && p[0] < 25) && (p[1] > 8 && p[1] < 27)) ||
            (p[0] > 25 && p[0] < 70) ||
            ((p[0] > -12 && p[0] < -1) && p[1] > 27));
  }
};

/**
 * Derived class for TPD-43 inclusions in the brain.
 */
template <int dim>
class TPD43_inclusions : public SeedingRegion<dim>
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

/**
 * Specialization of the TPD43_inclusions class for 3D problems.
 */
template <>
class TPD43_inclusions<3> : public SeedingRegion<3>
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

/**
 * Specialization of the TPD43_inclusions class for 2D problems.
 */
template <>
class TPD43_inclusions<2> : public SeedingRegion<2>
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

/**
 * Factory function to create seeding regions.
 * The seeding region is specified by the user in the main.cpp file.
 *
 * @param region Name of the seeding region.
 * @return Unique pointer to the seeding region.
 * @throw std::invalid_argument if the seeding region is invalid.
 */
template <int dim>
std::unique_ptr<SeedingRegion<dim>>
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

/**
 * Class to define the grey matter region in the brain.
 *
 * @tparam dim Dimension of the problem.
 */
template <int dim>
class Grey_matter
{
public:
  Grey_matter()
  {}

  /**
   * Check if a point is inside the grey matter region.
   *
   * @param p Point to check.
   * @return True if the point is inside the grey matter region, false otherwise.
   */
  static bool
  check_region(const Point<dim> &p)
  {
    return false;
  }
};

/**
 * Specialization of the Grey_matter class for 3D problems.
 */
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

/**
 * Specialization of the Grey_matter class for 2D problems.
 */
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
  // semi-axes length of the ellipse representing the grey matter region.
  static constexpr double a = 60.0;
  static constexpr double b = 40.0;
};

/**
 * Class to define the brain region where fibers are oriented in the
 * circumferential direction.
 *
 * @tparam dim Dimension of the problem.
 */
template <int dim>
class Axonal_region
{
public:
  Axonal_region()
  {}

  /**
   * Check if a point is inside the axonal region.
   *
   * @param p Point to check.
   * @return True if the point is inside the axonal region, false otherwise.
   */
  static bool
  check_region(const Point<dim> &p)
  {
    return false;
  }

  /**
   * Compute the radial direction at a given point.
   *
   * @param p Point at which to compute the radial direction.
   * @param global_center Global center of the brain.
   * @return Tensor<1, dim> representing the radial direction at the given point.
   */
  static Tensor<1, dim>
  compute_radial_direction(const Point<dim> &p, const Point<dim> &global_center)
  {
    return Tensor<1, dim>();
  }

  /**
   * Compute the circumferential direction at a given point.
   *
   * @param p Point at which to compute the circumferential direction.
   * @param global_center Global center of the brain.
   * @return Tensor<1, dim> representing the circumferential direction at the given point.
   */
  static Tensor<1, dim>
  compute_circumferential_direction(const Point<dim> &p,
                                    const Point<dim> &global_center)
  {
    return Tensor<1, dim>();
  }

  /**
   * Compute the axon-based direction at a given point.
   * The axon-based direction is the circumferential direction if the point is
   * inside the axonal region, and the radial direction otherwise.
   *
   * @param p Point at which to compute the axon-based direction.
   * @param global_center Global center of the brain.
   * @return Tensor<1, dim> representing the axon-based direction at the given point.
   */
  static Tensor<1, dim>
  compute_axon_based_direction(const Point<dim> &p,
                               const Point<dim> &global_center)
  {
    if (Axonal_region<dim>::check_region(p))
      {
        return compute_circumferential_direction(p, global_center);
      }
    else
      return compute_radial_direction(p, global_center);
  }

  /**
   * Get the axonal direction at a given point.
   * The axonal direction can be either radial, circumferential, or axon-based.
   *
   * @param p Point at which to compute the axonal direction.
   * @param global_center Global center of the brain.
   * @param orientation Orientation of the axonal direction.
   * @return Tensor<1, dim> representing the axonal direction at the given point.
   */
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

/**
 * Specialization of the Axonal_region class for 3D problems.
 */
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

    Tensor<1, 3> arbitrary_vector;
    if (std::abs(radial[0]) < std::abs(radial[1]) &&
        std::abs(radial[0]) < std::abs(radial[2]))
      {
        arbitrary_vector[0] = 1.0;
        arbitrary_vector[1] = 0.0;
        arbitrary_vector[2] = 0.0;
      }
    else if (std::abs(radial[1]) < std::abs(radial[2]))
      {
        arbitrary_vector[0] = 0.0;
        arbitrary_vector[1] = 1.0;
        arbitrary_vector[2] = 0.0;
      }
    else
      {
        arbitrary_vector[0] = 0.0;
        arbitrary_vector[1] = 0.0;
        arbitrary_vector[2] = 1.0;
      }

    // Compute the azimuthal direction as the cross product of radial and
    // arbitrary_vector
    Tensor<1, 3> azimuthal;
    azimuthal[0] =
      radial[1] * arbitrary_vector[2] - radial[2] * arbitrary_vector[1];
    azimuthal[1] =
      radial[2] * arbitrary_vector[0] - radial[0] * arbitrary_vector[2];
    azimuthal[2] =
      radial[0] * arbitrary_vector[1] - radial[1] * arbitrary_vector[0];

    azimuthal /= azimuthal.norm();

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

/**
 * Specialization of the Axonal_region class for 2D problems.
 */
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


template <int dim>
class RadialDirection : Function<dim>
{
public:
  RadialDirection(const Point<dim> &global_center_)
    : global_center(global_center_)
  {}

  virtual void
  vector_value(const Point<dim> &p, Vector<double> &values) const override
  {
    Tensor<1, dim> radial = p - global_center;
    radial /= radial.l2_norm();
    values[0] = radial[0];
    values[1] = radial[1];
    if (dim == 3)
      values[2] = radial[2];
  }

  virtual double
  value(const Point<dim> &p, const unsigned int component = 0) const override
  {
    Tensor<1, dim> radial = p - global_center;
    radial /= radial.l2_norm();
    return radial[component];
  }

protected:
  const Point<dim> global_center;
};

template <int dim>
class CirtumferentialDirection : Function<dim>
{
public:
  CirtumferentialDirection(const Point<dim> &global_center_)
    : global_center(global_center_)
  {}

  virtual void
  vector_value(const Point<dim> &p, Vector<double> &values) const override
  {
    Vector<double> radial = p - global_center;
    if constexpr (dim == 2)
      {
        Vector<double> azimuthal;
        azimuthal[0] = -p[1];
        azimuthal[1] = p[0];
        azimuthal /= azimuthal.l2_norm();

        // Cross product between radial and azimuthal directions
        values[0] = radial[1] * azimuthal[0] - radial[0] * azimuthal[1];
        values[1] = radial[0] * azimuthal[0] + radial[1] * azimuthal[1];
      }
    else
      {
        Vector<double> arbitrary_vector;
        if (std::abs(radial[0]) < std::abs(radial[1]) &&
            std::abs(radial[0]) < std::abs(radial[2]))
          {
            arbitrary_vector[0] = 1.0;
            arbitrary_vector[1] = 0.0;
            arbitrary_vector[2] = 0.0;
          }
        else if (std::abs(radial[1]) < std::abs(radial[2]))
          {
            arbitrary_vector[0] = 0.0;
            arbitrary_vector[1] = 1.0;
            arbitrary_vector[2] = 0.0;
          }
        else
          {
            arbitrary_vector[0] = 0.0;
            arbitrary_vector[1] = 0.0;
            arbitrary_vector[2] = 1.0;
          }

        // Compute the azimuthal direction as the cross product of radial and
        // arbitrary_vector
        Vector<double> azimuthal;
        azimuthal[0] =
          radial[1] * arbitrary_vector[2] - radial[2] * arbitrary_vector[1];
        azimuthal[1] =
          radial[2] * arbitrary_vector[0] - radial[0] * arbitrary_vector[2];
        azimuthal[2] =
          radial[0] * arbitrary_vector[1] - radial[1] * arbitrary_vector[0];

        azimuthal /= azimuthal.l2_norm();

        values[0] = radial[1] * azimuthal[2] - radial[2] * azimuthal[1];
        values[1] = radial[2] * azimuthal[0] - radial[0] * azimuthal[2];
        values[2] = radial[0] * azimuthal[1] - radial[1] * azimuthal[0];
      }
  }

  virtual double
  value(const Point<dim> &p, const unsigned int component = 0) const override
  {
    Vector<double> radial = p - global_center;
    Vector<double> circumferential;
    if constexpr (dim == 2)
      {
        Vector<double> azimuthal;
        azimuthal[0] = -p[1];
        azimuthal[1] = p[0];
        azimuthal /= azimuthal.l2_norm();

        // Cross product between radial and azimuthal directions
        circumferential[0] =
          radial[1] * azimuthal[0] - radial[0] * azimuthal[1];
        circumferential[1] =
          radial[0] * azimuthal[0] + radial[1] * azimuthal[1];
      }
    else
      {
        Vector<double> arbitrary_vector;
        if (std::abs(radial[0]) < std::abs(radial[1]) &&
            std::abs(radial[0]) < std::abs(radial[2]))
          {
            arbitrary_vector[0] = 1.0;
            arbitrary_vector[1] = 0.0;
            arbitrary_vector[2] = 0.0;
          }
        else if (std::abs(radial[1]) < std::abs(radial[2]))
          {
            arbitrary_vector[0] = 0.0;
            arbitrary_vector[1] = 1.0;
            arbitrary_vector[2] = 0.0;
          }
        else
          {
            arbitrary_vector[0] = 0.0;
            arbitrary_vector[1] = 0.0;
            arbitrary_vector[2] = 1.0;
          }

        // Compute the azimuthal direction as the cross product of radial and
        // arbitrary_vector
        Vector<double> azimuthal;
        azimuthal[0] =
          radial[1] * arbitrary_vector[2] - radial[2] * arbitrary_vector[1];
        azimuthal[1] =
          radial[2] * arbitrary_vector[0] - radial[0] * arbitrary_vector[2];
        azimuthal[2] =
          radial[0] * arbitrary_vector[1] - radial[1] * arbitrary_vector[0];

        circumferential[0] =
          radial[1] * azimuthal[2] - radial[2] * azimuthal[1];
        circumferential[1] =
          radial[2] * azimuthal[0] - radial[0] * azimuthal[2];
        circumferential[2] =
          radial[0] * azimuthal[1] - radial[1] * azimuthal[0];
      }
    return circumferential[component];
  }

protected:
  const Point<dim> global_center;
};


/*
class AxonBasedDirection : Function<dim>
{
public:
  AxonBasedDirection(const Point<dim>  &global_center_,
                     const std::string &orientation_)
    : global_center(global_center_)
    , orientation(orientation_)
  {}

  virtual void
  vector_value(const Point<dim> &p, Vector<double> &values) const override
  {
    if (Axonal_region<dim>::check_region(p))
      {
        CirtumferentialDirection<dim> circumferential(global_center);
        circumferential.vector_value(p, values);
      }
    else
      {
        RadialDirection<dim> radial(global_center);
        radial.vector_value(p, values);
      }
  }

  virtual double
  value(const Point<dim> &p, const unsigned int component = 0) const override
  {
    if (Axonal_region<dim>::check_region(p))
      {
        CirtumferentialDirection<dim> circumferential(global_center);
        return circumferential.value(p, component);
      }
    else
      {
        RadialDirection<dim> radial(global_center);
        return radial.value(p, component);
      }
  }

protected:
  const Point<dim>  global_center;
  const std::string orientation;
};
*/
#endif