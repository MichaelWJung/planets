#pragma once

#include "simulator.h"

#include <cstddef>
#include <deque>
#include <vector>

namespace planets {

class Trails {
public:
  explicit Trails(std::size_t capacity);

  void record(const std::vector<Body> &bodies);

  [[nodiscard]] const std::vector<std::deque<position_t>> &data() const;

private:
  std::size_t capacity_;
  std::vector<std::deque<position_t>> trails_;
};

} // namespace planets
