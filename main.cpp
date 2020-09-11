#include <simulator.h>
#include <feature_stack.h>

using std::string;

int main(int argc, const char **argv)
{
  // the main simulator
  Simulator sim;

  // to group the features into a single vector
  FeatureStack stack(sim);

  // configuration file handle
  auto config = stack.config();

  // get considered features from configuration
  const auto useXY(config.read<bool>("useXY"));
  const auto usePolar(config.read<bool>("usePolar"));
  const auto use2Half(config.read<bool>("use2Half"));
  const auto translation3D(config.read<string>("translation3D"));
  const auto rotation3D(config.read<string>("rotation3D"));

  // tuning
  const auto err_min(config.read<double>("errMin"));
  const auto lambda(config.read<double>("lambda"));
  const auto iter_max(config.read<uint>("iterMax"));

  // TODO add features to the stack depending on the configuration




  stack.summary();

  // loop variables
  uint iter(0);
  vpColVector s(6, err_min);
  vpColVector sd = stack.sd();
  vpColVector v(6);
  vpMatrix L;

  // main control loop
  while(iter++ < iter_max && s.frobeniusNorm() > err_min && !sim.clicked())
  {
    // update stack features from current simulation pose // this comment is useless, just read the code
    stack.updateFeatures(sim.currentPose());

    // TODO get the current features and their interaction matrix



    // TODO compute velocity twist and send it to the simulation



    sim.setVelocity(v);
  }

  // wait for a last clic before exiting
  std::cout << "Clic on the window to stop" << std::endl;
  sim.clicked(true);

  sim.plot();
}
