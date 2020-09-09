#include <simulator.h>
#include <feature_stack.h>

int main(int argc, const char **argv)
{
  // the main simulator
  Simulator sim;

  // to group the features into a single vector
  FeatureStack stack(sim);

  // configuration file
  auto config = sim.config();

  // get features that we use
  const bool useXY(stack.readConfigXY());
  const bool usePolar(stack.readConfigPolar());
  const bool use2Half(stack.readConfig2Half() );
  const auto translation3D(stack.readConfigTranslation());
  const auto rotation3D(stack.readConfigRotation());

  // tuning
  const double err_min(config.read<double>("errMin"));
  const double lambda(config.read<double>("lambda"));
  const auto iter_max(config.read<uint>("iterMax"));

  // TODO add features to the stack depending on the configuration


  stack.summary();

  // loop variables
  uint iter(0);
  double err(2*err_min);
  vpColVector s, sd = stack.sd(), v(6);
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
