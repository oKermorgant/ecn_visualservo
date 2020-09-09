#include <simulator.h>
#include <visp/vpMeterPixelConversion.h>

template <typename... Args> inline void UNUSED(Args&&...) {}

Simulator::Simulator() : config_manager(std::string(BASE_PATH) + "config.yaml")
{
  UNUSED(system("killall python"));

  // init images and display
  Iint.init(480, 640, 255);
  Iext.init(480, 640, 255);
  Dint.init(Iint, 100, 100, "Internal view");
  Dext.init(Iext, 100, 100, "External view");
  vpDisplay::setWindowPosition(Iint, 0, 0);
  vpDisplay::setWindowPosition(Iext, 700, 0);
  vpDisplay::display(Iint);
  vpDisplay::flush(Iint);
  vpDisplay::display(Iext);
  vpDisplay::flush(Iext);

  // config file
  // initial pose
  const auto start = config_manager.read<std::string>("startPos");
  const auto end = config_manager.read<std::string>("endPos");
  cMo = config_manager.read<vpHomogeneousMatrix>(start);
  cdMo = config_manager.read<vpHomogeneousMatrix>(end);
  robot.setPosition(cMo);
  robot.setMaxTranslationVelocity(10);
  robot.setMaxRotationVelocity(10);

  // simulator
  sim.initScene ( vpWireFrameSimulator::SQUARE_10CM, vpWireFrameSimulator::D_STANDARD );
  cam.initPersProjWithoutDistortion(1000, 1000, 320, 240);
  sim.setCurrentViewColor(vpColor::blue);
  sim.setDesiredViewColor(vpColor::red);
  sim.setCameraColor(vpColor::green);
  sim.setCameraPositionRelObj(cMo);
  sim.setDesiredCameraPosition(cdMo);
  sim.setExternalCameraPosition(vpHomogeneousMatrix(0.0, 0, 4.5, vpMath::rad(0), vpMath::rad(-30), 0));
  sim.setInternalCameraParameters(cam);
  sim.setExternalCameraParameters(cam);

  // 3D points
  const double ps = 0.05;
  points[0].setWorldCoordinates ( -ps,-ps,0 );
  points[3].setWorldCoordinates ( -ps,ps,0 );
  points[2].setWorldCoordinates ( ps,ps,0 );
  points[1].setWorldCoordinates ( ps,-ps,0 );

  // CoG
  double x(0), y(0), z(0);
  for(const auto &p: points)
  {
    x += p.get_oX();
    y += p.get_oY();
    z += p.get_oZ();
  }
  center = vpPoint(x/4, y/4, z/4);
}

void Simulator::setVelocity(const vpColVector &v)
{
  robot.setVelocity(vpRobotCamera::CAMERA_FRAME, v);
  vel = v;
  robot.getPosition(cMo);
  sim.setCameraPositionRelObj ( cMo ) ;
  sim.getInternalImage ( Iint );

  pose.buildFrom(cMo);

  uint idx(0);
  for(auto &P: points)
  {
    P.track(cMo);
    uv[2*idx] = computeU(P);
    uv[2*idx+1] = computeV(P);
    idx++;
    vpImagePoint ip;
    vpMeterPixelConversion::convertPoint(cam, P.get_x(), P.get_y(), ip);
    history.push_back(ip);
  }
  center.track(cMo);
  uv[2*idx] = computeU(center);
  uv[2*idx+1] = computeV(center);
  vpImagePoint ip;
  vpMeterPixelConversion::convertPoint(cam, center.get_x(), center.get_y(), ip);
  history.push_back(ip);

  for(const auto &ip: history)
    vpDisplay::displayPoint(Iint, ip, vpColor::darkRed, 2);
  vpDisplay::flush ( Iint );
  vpDisplay::display ( Iint );
  sim.getExternalImage(Iext);
  vpDisplay::flush ( Iext );
  vpDisplay::display ( Iext );
  vpTime::wait(t0, dt);
  t0 = vpTime::measureTimeSecond();

  logger.update();
}

void Simulator::initLog(const std::string &exp_id, const std::string &legend)
{
  std::string log_dir(BASE_PATH);
  log_dir += "results/" + exp_id + "/";
  config_manager.setDirName(log_dir);
  logger.setSavePath(log_dir);

  t0 = vpTime::measureTimeSecond();
  // 2D points XY (+center)
  // image limits
  const double w(Iint.getCols());
  const double h(Iint.getRows());

  uv.resize(10);
  logger.saveXY(uv, "image", "[P_1, P_2, P_3, P_4, CoG]", "u", "v");
  logger.setLineType("[C0,C1,C2,C3,C4]");
  logger.showFixedObject({{0,0},{0,h},{w,h},{w,0}}, "[[0,1],[1,2],[2,3],[3,0]]", "k-");

  std::vector<std::vector<double>> des, cur;
  for(auto &P: points)
  {
    P.track(cdMo);
    des.push_back({computeU(P), computeV(P)});
    P.track(cMo);
    cur.push_back({computeU(P), computeV(P)});
  }
  logger.showFixedObject(des, "[[0,1],[1,2],[2,3],[3,0]]", "rs--");
  logger.showFixedObject(cur, "[[0,1],[1,2],[2,3],[3,0]]", "bs--");
  logger.setPlotArgs("--xLim -10 " + std::to_string(w+10) + " --yLim -10 " + std::to_string(h+10));

  // 3D pose
  logger.save3Dpose(pose, "pose", "'"+legend+"'", true);
  vpPoseVector pose_d(cdMo);
  logger.showMovingCamera({pose_d[0], pose_d[1], pose_d[2], pose_d[3], pose_d[4], pose_d[5]});
  logger.showFixedRectangle(-.05, -.05, .05, .05, "b");

  // velocity
  vel.resize(6);
  logger.save(vel, "v", "[v_x,v_y,v_z,\\omega_x,\\omega_y,\\omega_z]", "Velocity");
}
