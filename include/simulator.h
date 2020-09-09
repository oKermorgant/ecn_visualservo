#ifndef VSSIMULATOR_H
#define VSSIMULATOR_H

#include <visp3/core/vpCameraParameters.h>
#include <visp3/core/vpHomogeneousMatrix.h>
#include <visp3/core/vpImage.h>
#include <visp3/core/vpMath.h>
#include <visp3/gui/vpDisplayX.h>
#include <visp3/io/vpImageIo.h>
#include <visp/vpRobotCamera.h>
#include <visp3/robot/vpWireFrameSimulator.h>
#include <log2plot/config_manager.h>
#include <log2plot/logger.h>
#include <chrono>

class FeatureStack;

class Simulator
{
public:
  Simulator();
  log2plot::ConfigManager config() const
  {
    return config_manager;
  }
  void setVelocity(const vpColVector &v);

  vpHomogeneousMatrix currentPose() const
  {
    return cMo;
  }
  vpHomogeneousMatrix desiredPose() const
  {
    return cdMo;
  }
  vpPoint cog() const
  {
    return center;
  }
  std::array<vpPoint, 4> observedPoints() const
  {
    return points;
  }

  void plot()
  {
    logger.plot(true);
    config_manager.saveConfig();
  }

  bool clicked(bool wait = false)
  {
    return  vpDisplay::getClick(Iint,wait);
  }

protected:
  friend class FeatureStack;
  vpImage<vpRGBa> Iint, Iext;
  vpDisplayX Dint, Dext;
  vpWireFrameSimulator sim;
  vpCameraParameters cam;
  vpRobotCamera robot;
  std::vector<vpImagePoint> history;

  std::array<vpPoint, 4> points;
  vpPoint center;
  vpColVector uv, vel;

  // time sampling
  const double dt = 0.01 * 1000;  // ms
  double t0;


  // from configuration
  vpHomogeneousMatrix cMo, cdMo;
  vpPoseVector pose;
  log2plot::ConfigManager config_manager;
  log2plot::Logger logger;

  void initLog(const std::string &exp_id, const std::string &legend);
  double computeV(const vpPoint &P) const
  {
    return Iint.getRows() - (cam.get_v0() + cam.get_py()*P.get_y());
  }
  double computeU(const vpPoint &P) const
  {
    return cam.get_u0() + cam.get_px()*P.get_x();
  }
};

#endif // VSSIMULATOR_H
