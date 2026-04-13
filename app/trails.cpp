#include "trails.h"

namespace planets {

Trails::Trails(std::size_t capacity) : capacity_(capacity) {}

void Trails::record(const std::vector<Body> &bodies) {
  if (trails_.size() != bodies.size())
    trails_.resize(bodies.size());

  for (std::size_t i = 0; i < bodies.size(); ++i) {
    trails_[i].push_back(bodies[i].position);
    if (trails_[i].size() > capacity_)
      trails_[i].pop_front();
  }
}

const std::vector<std::deque<position_t>> &Trails::data() const {
  return trails_;
}

} // namespace planets
