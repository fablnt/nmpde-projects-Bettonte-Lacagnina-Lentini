#include <memory>
#include <string>

// SeedingRegion is a base struct for different seeding regions
struct SeedingRegion
{
public:
  SeedingRegion(const double x_min_,
                const double x_max_,
                const double y_min_,
                const double y_max_,
                const double z_min_,
                const double z_max_)
    : x_min(x_min_)
    , x_max(x_max_)
    , y_min(y_min_)
    , y_max(y_max_)
    , z_min(z_min_)
    , z_max(z_max_)
  {}

  virtual ~SeedingRegion()
  {}

  const double x_min;
  const double x_max;
  const double y_min;
  const double y_max;
  const double z_min;
  const double z_max;
};
;

// todo
struct Tau_inclusions : public SeedingRegion
{
public:
  Tau_inclusions()
    : SeedingRegion(63.0, 81.0, 60.0, 90.0, 46.0, 67.0)
  {}
};

// todo
struct Amyloid_Beta_deposits : public SeedingRegion
{
public:
  Amyloid_Beta_deposits()
    : SeedingRegion(43.0, 80.0, 22.0, 90.0, 80.0, 115.0)
  {}
};

// todo
struct TPD43_inclusions : public SeedingRegion
{
public:
  TPD43_inclusions()
    : SeedingRegion(52.0, 80.0, 73.0, 113.0, 41.0, 115.0)
  {}
};

// Factory function to create seeding regions.
inline std::unique_ptr<SeedingRegion>
getSeedingRegion(std::string region)
{
  if (region == "Tau inclusions")
    return std::make_unique<Tau_inclusions>();
  else if (region == "Amyloid-Beta deposits")
    return std::make_unique<Amyloid_Beta_deposits>();
  else if (region == "TPD-43 inclusions")
    return std::make_unique<TPD43_inclusions>();
  else
    throw std::invalid_argument("Invalid seeding region");
}
