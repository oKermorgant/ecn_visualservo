#include <simulator.h>
#include <feature_stack.h>

int main(int argc, const char **argv)
{
  // the main simulator
  Simulator sim;

  // to group the features into a single vector
  FeatureStack stack(sim);

  // configuration file
  auto config = stack.config();

  // get considered features from configuration
  const auto useXY(stack.readConfigXY());
  const auto usePolar(stack.readConfigPolar());
  const auto use2Half(stack.readConfig2Half() );
  const auto translation3D(stack.readConfigTranslation());
  const auto rotation3D(stack.readConfigRotation());

  // tuning
  const auto err_min(config.read<double>("errMin"));
  const auto lambda(config.read<double>("lambda"));
  const auto iter_max(config.read<uint>("iterMax"));

  // TODO add features to the stack depending on the configuration




  stack.summary();

  // loop variables
  uint iter(0);
  auto err(2*err_min);
  vpColVector s(6, 1), sd = stack.sd(), v(6);
  vpMatrix L;
  vpHomogeneousMatrix cMo;

  while(iter++ < iter_max && err > err_min && !sim.clicked())
  {
    // current transform
    cMo = sim.currentPose();
    stack.updateFeatures(cMo);

    // TODO get the current features and their interaction matrix



    // TODO compute velocity twist and send it to the simulation



    sim.setVelocity(v);

    // register this error
    err = s.frobeniusNorm();
  }
  if(iter == iter_max || err < err_min)
  {
    std::cout << "Clic on the window to stop" << std::endl;
    sim.clicked(true);
  }

  sim.plot();
}
