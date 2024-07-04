#ifndef REGIONS_HPP
#define REGIONS_HPP

#include <memory>
#include <string>

using namespace dealii;
/**
 * Base class for initial seeding regions in the brain.
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
 * @throw std::invalid_argument if the seeding region name is invalid.
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
  bool
  check_grey_matter(const Point<dim> &p)
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

  bool
  check_grey_matter(const Point<3> &p)
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

  bool
  check_grey_matter(const Point<2> &p)
  {
    return (p[0] * p[0]) / (a * a) + (p[1] * p[1]) / (b * b) > 1.0;
  }

private:
  // semi-axes length of the ellipse representing the grey matter region in 2D.
  static constexpr double a = 60.0;
  static constexpr double b = 40.0;
};


/**
 * Class to define the Radial direction of the fibers in the brain, with respect
 * to the global center of the domain.
 *
 * @tparam dim Dimension of the problem.
 */
template <int dim>
class RadialDirection : public Function<dim>
{
public:
  RadialDirection(const Point<dim> &global_center_)
    : global_center(global_center_)
  {}

  virtual void
  vector_value(const Point<dim> &p, Vector<double> &values) const override
  {
    Vector<double> radial(dim);
    for (unsigned int i = 0; i < dim; ++i)
      radial[i] = p[i] - global_center[i];
    radial /= radial.l2_norm();

    values[0] = radial[0];
    values[1] = radial[1];
    if (dim == 3)
      values[2] = radial[2];
  }

  virtual double
  value(const Point<dim> &p, const unsigned int component = 0) const override
  {
    Vector<double> radial(dim);
    for (unsigned int i = 0; i < dim; ++i)
      radial[i] = p[i] - global_center[i];
    radial /= radial.l2_norm();

    return radial[component];
  }

protected:
  const Point<dim> global_center;
};
/**
 * Class to define the Circumferential direction of the fibers in the brain
 * relative to the global center of the domain. The cirumferential direction is
 * tangent to the radial direction.
 *
 * In the 2D problem, the circumferential direction is obtained by rotating
 * the radial vector by 90° clockwise.
 *
 * In the 3D problem, the circumferential direction is defined using cross
 * products:
 * First, an orthogonal vector is computed as the cross product between the
 * radial vector and an arbitrary vector. The arbitrary vector is chosen to be
 * as orthogonal as possible to the radial direction.
 * Second, the circumferential direction is computed as the cross product
 * between the radial vector and the orthogonal vector. This ensures that the
 * circumferential direction is tangent to the radial direction.
 *
 * @tparam dim Dimension of the problem.
 */
template <int dim>
class CircumferentialDirection : public Function<dim>
{
public:
  CircumferentialDirection(const Point<dim> &global_center_)
    : global_center(global_center_)
  {}

  virtual void
  vector_value(const Point<dim> &p, Vector<double> &values) const override
  {
    Vector<double> radial(dim);
    for (unsigned int i = 0; i < dim; ++i)
      radial[i] = p[i] - global_center[i];
    radial /= radial.l2_norm();
    Vector<double> orthogonal(dim);
    if constexpr (dim == 2)
      {
        orthogonal[0] = -p[1];
        orthogonal[1] = p[0];
        orthogonal /= orthogonal.l2_norm();
      }
    else
      {
        Vector<double> arbitrary_vector(dim);
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

        // Compute the orthogonal direction as the cross product of radial and
        // arbitrary_vector.
        orthogonal[0] =
          radial[1] * arbitrary_vector[2] - radial[2] * arbitrary_vector[1];
        orthogonal[1] =
          radial[2] * arbitrary_vector[0] - radial[0] * arbitrary_vector[2];
        orthogonal[2] =
          radial[0] * arbitrary_vector[1] - radial[1] * arbitrary_vector[0];

        orthogonal /= orthogonal.l2_norm();

        // Compute the circumferential direction as the cross product of radial
        // and orthogonal directions.
        values[0] = radial[1] * orthogonal[2] - radial[2] * orthogonal[1];
        values[1] = radial[2] * orthogonal[0] - radial[0] * orthogonal[2];
        values[2] = radial[0] * orthogonal[1] - radial[1] * orthogonal[0];
      }
  }

  virtual double
  value(const Point<dim> &p, const unsigned int component = 0) const override
  {
    Vector<double> radial(dim);
    for (unsigned int i = 0; i < dim; ++i)
      radial[i] = p[i] - global_center[i];
    radial /= radial.l2_norm();
    Vector<double> circumferential(dim);
    Vector<double> orthogonal(dim);
    if constexpr (dim == 2)
      {
        orthogonal[0] = -p[1];
        orthogonal[1] = p[0];
        orthogonal /= orthogonal.l2_norm();
      }
    else
      {
        Vector<double> arbitrary_vector(dim);
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

        // Compute the orthogonal direction as the cross product of radial and
        // arbitrary_vector.
        Vector<double> orthogonal;
        orthogonal[0] =
          radial[1] * arbitrary_vector[2] - radial[2] * arbitrary_vector[1];
        orthogonal[1] =
          radial[2] * arbitrary_vector[0] - radial[0] * arbitrary_vector[2];
        orthogonal[2] =
          radial[0] * arbitrary_vector[1] - radial[1] * arbitrary_vector[0];

        orthogonal /= orthogonal.l2_norm();

        // Compute the circumferential direction as the cross product of radial
        // and orthogonal directions.
        circumferential[0] =
          radial[1] * orthogonal[2] - radial[2] * orthogonal[1];
        circumferential[1] =
          radial[2] * orthogonal[0] - radial[0] * orthogonal[2];
        circumferential[2] =
          radial[0] * orthogonal[1] - radial[1] * orthogonal[0];
      }
    return circumferential[component];
  }

protected:
  const Point<dim> global_center;
};

/**
 * Class to define the direction of the fibers in the brain following those of
 * the axons. Based on the locations of the points, the direction is either
 * radial or circumferential.
 *
 * @tparam dim Dimension of the problem.
 */
template <int dim>
class AxonBasedDirection : public Function<dim>
{
public:
  AxonBasedDirection(const Point<dim> &global_center_)
    : global_center(global_center_)
    , circumferential(global_center_)
    , radial(global_center_)
  {}

  virtual void
  vector_value(const Point<dim> &p, Vector<double> &values) const override
  {
    if (check_axonal_region(p))
      circumferential.vector_value(p, values);
    else
      radial.vector_value(p, values);
  }

  virtual double
  value(const Point<dim> &p, const unsigned int component = 0) const override
  {
    if (check_axonal_region(p))
      return circumferential.value(p, component);
    else
      return radial.value(p, component);
  }

  /**
   * Check in which region the point is located. If the points is in the region
   * where the axonal fiber orientation is circumferential, return true.
   * Otherwise, the point is in the region where the axonal fiber orientation is
   * radial and return false.
   *
   * @param p Point to check.
   * @return True if the point is in the region where the axonal fiber orientation
   *         is circumferential, false otherwise.
   */
  bool
  check_axonal_region(const Point<dim> &p) const
  {
    if constexpr (dim == 3)
      {
        return (p[0] < 60 && p[0] > 40) && (p[1] < 110 && p[1] > 34) &&
               (p[2] < 80 && p[2] > 50);
      }
    else if constexpr (dim == 2)
      {
        double a = 30.0;
        double b = 20.0;
        return (p[0] * p[0]) / (a * a) + (p[1] * p[1]) / (b * b) < 1.0;
      }
    else
      {
        static_assert(dim == 2 || dim == 3, "Unsupported dimension");
        return false; // Should not reach here
      }
  }

protected:
  const Point<dim>              global_center;
  CircumferentialDirection<dim> circumferential;
  RadialDirection<dim>          radial;
};

/**
 * Factory function to create the direction of the fibers in the brain.
 * The direction is specified by the user in the main.cpp file.
 *
 * @param orientation Name of the orientation.
 * @param global_center Global center of the brain.
 * @return Unique pointer to the direction.
 * @throw std::invalid_argument if the orientation name is invalid.
 */
template <int dim>
std::unique_ptr<Function<dim>>
get_direction(std::string orientation, const Point<dim> &global_center)
{
  if (orientation == "radial")
    return std::make_unique<RadialDirection<dim>>(global_center);
  else if (orientation == "circumferential")
    return std::make_unique<CircumferentialDirection<dim>>(global_center);
  else if (orientation == "axon-based")
    return std::make_unique<AxonBasedDirection<dim>>(global_center);
  else
    throw std::invalid_argument("Invalid orientation");
}

#endif