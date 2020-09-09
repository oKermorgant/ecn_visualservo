#include <feature_stack.h>
#include <logger.h>

log2plot::Logger Logger::logger(std::string(BASE_PATH));

void recap(const std::string s, int n, int dim)
{
  if(n)
    std::cout << s << ": " << n << " (dim " << dim*n << ")" << std::endl;
}

void FeatureStack::summary() const
{
  std::cout << "Total feature dimension: " << s_rows << std::endl;
  recap(" - Points XY", pointsXY.size(), 2);
  recap(" - Points Polar", pointsPolar.size(), 2);
  recap(" - Points depths", depths.size(), 1);
  if(translation.des != TranslationDescriptor::NONE)
    recap(" - 3D translations", 1, 3);
  if(rotation.des != RotationDescriptor::NONE)
    recap(" - 3D rotations", 1, 3);
}

bool FeatureStack::readConfigXY() const
{
  return simulator.config_manager.read<bool>("useXY");
}
bool FeatureStack::readConfigPolar() const
{
  return simulator.config_manager.read<bool>("usePolar");
}
bool FeatureStack::readConfig2Half() const
{
  return simulator.config_manager.read<bool>("use2Half");
}
TranslationDescriptor FeatureStack::readConfigTranslation() const
{
  const auto des(simulator.config_manager.read<std::string>("translation3D"));
  if(des == "cTo")
    return TranslationDescriptor::cTo;
  else if(des == "cdTc")
    return TranslationDescriptor::cdTc;
  return TranslationDescriptor::NONE;
}
RotationDescriptor FeatureStack::readConfigRotation() const
{
  const auto des(simulator.config_manager.read<std::string>("rotation3D"));
  if(des == "cdRc")
    return RotationDescriptor::cdRc;
  else if(des == "cRcd")
    return RotationDescriptor::cRcd;
  return RotationDescriptor::NONE;
}

vpColVector FeatureStack::sd()
{
  if(sd_.size() == 0)
  {
    sd_.resize(s_rows);
    computeFeatures(cdMo, false);
  }
  return sd_;
}

void FeatureStack::updateFeatures(const vpHomogeneousMatrix &cMo)
{
  if(!init_done)
    initLog();
  computeFeatures(cMo, true);
  e_ = s_-sd();
}

void FeatureStack::computeFeatures(const vpHomogeneousMatrix &cMo, bool current)
{
  uint row(0);
  if(translation.des != TranslationDescriptor::NONE)
  {
    if(translation.des == TranslationDescriptor::cTo)
      translation.feature.buildFrom(cMo);
    else
      translation.feature.buildFrom(cdMo*cMo.inverse());
    //std::cout << "T = " << translation.feature.get_s().t() << std::endl;
    //std::cout << "Lt = \n" << translation.feature.interaction() << std::endl;
    row = update(row, translation.feature, current);
  }

  if(rotation.des != RotationDescriptor::NONE)
  {
    if(rotation.des == RotationDescriptor::cRcd)
      rotation.feature.buildFrom(cMo*cdMo.inverse());
    else
      rotation.feature.buildFrom(cdMo*cMo.inverse());
    row = update(row, rotation.feature, current);
  }

  auto xy(pointsXY.begin());
  auto rt(pointsPolar.begin());
  auto d(depths.begin());
  for(auto &[P, descriptor, zd]: points3D)
  {
    P.track(cMo);
    const double z(P.get_Z());
    // Z-estimation
    if(current)
    {
      if(z_estim > 0)
        P.set_Z(z_estim);
      else if(z_estim == 0)
        P.set_Z(zd);
    }
    if(descriptor == PointDescriptor::XY)
    {
      xy->buildFrom(P.get_x(), P.get_y(), P.get_Z());
      row = update(row, *xy, current);
      xy++;
    }
    else if(descriptor == PointDescriptor::Polar)
    {
      vpFeatureBuilder::create(*rt, P);
      row = update(row, *rt, current);
      rt++;
    }
    else
    {
      d->buildFrom(P.get_x(), P.get_y(), P.get_Z(), log(z/zd));
      row = update(row, *d, current);
      d++;
    }
  }
}

void FeatureStack::addFeaturePoint(vpPoint P, PointDescriptor descriptor)
{
  if(descriptor == PointDescriptor::XY)
  {
    pointsXY.push_back(vpFeaturePoint());
    s_rows += 2;
  }
  else if(descriptor == PointDescriptor::Polar)
  {
    pointsPolar.push_back(vpFeaturePointPolar());
    s_rows += 2;
  }
  else
  {
    depths.push_back(vpFeatureDepth());
    s_rows += 1;
  }
  P.track(cdMo);
  points3D.push_back({P, descriptor, P.get_Z()});
}

void FeatureStack::setTranslation3D(TranslationDescriptor descriptor)
{
  if(descriptor == TranslationDescriptor::NONE)
    return;

  if(descriptor == TranslationDescriptor::cTo)
    translation.feature.setFeatureTranslationType(vpFeatureTranslation::cMo);
  else
    translation.feature.setFeatureTranslationType(vpFeatureTranslation::cdMc);

  s_rows += 3;
  translation.des = descriptor;
}

void FeatureStack::setRotation3D(RotationDescriptor descriptor)
{
  if(descriptor == RotationDescriptor::NONE)
    return;

  if(descriptor == RotationDescriptor::cRcd)
    rotation.feature.setFeatureThetaURotationType(vpFeatureThetaU::cRcd);
  else
    rotation.feature.setFeatureThetaURotationType(vpFeatureThetaU::cdRc);

  s_rows += 3;
  rotation.des = descriptor;
}

void FeatureStack::initLog()
{
  init_done = true;

  // build sub-dir and legend
  std::stringstream ss;
  std::string exp_id, legend3D;
  auto updatePath([&]()
  {
    const auto key(ss.str());
    auto legend_key(key);
    /* if(legend_key == "cdRc")
      legend_key = "{}^{c*}\\mathbf{R}_c";
    else if(legend_key == "cRcd")
      legend_key = "{}^c\\mathbf{R}_{c*}";
    else if(legend_key == "cTo")
        legend_key = "{}^c\\mathbf{T}_o";
    else if(legend_key == "cdTc")
      legend_key = "{}^{c*}\\mathbf{T}_c";*/
    if(exp_id == "")
    {
      exp_id = key;
      legend3D = key;
    }
    else
    {
      exp_id += "_" + key;
      legend3D += "+" + legend_key;
    }
    ss.str("");
  });

  std::stringstream legend;
  legend << "[";
  if(pointsXY.size())
  {
    ss << "XY" << pointsXY.size();
    updatePath();
    for(size_t i = 0; i < pointsXY.size(); ++i)
      legend << "'x_{" << i+1 << "}', 'y_{" << i+1 << "}',";
  }
  if(pointsPolar.size())
  {
    ss << "Polar" << pointsPolar.size();
    updatePath();
    for(size_t i = 0; i < pointsPolar.size(); ++i)
      legend << "'\\rho_{" << i+1 << "}', '\\theta_{" << i+1 << "}',";
  }
  if(depths.size())
  {
    ss << "Depth" << depths.size();
    updatePath();
    for(size_t i = 0; i < depths.size(); ++i)
      legend << "'Z_{" << i+1 << "}',";
  }
  if(translation.des != TranslationDescriptor::NONE)
  {
    ss << (translation.des == TranslationDescriptor::cTo ? "cTo" : "cdTc");
    updatePath();
    legend << "t_x, t_y, t_z,";
  }
  if(rotation.des != RotationDescriptor::NONE)
  {
    ss << (rotation.des == RotationDescriptor::cRcd ? "cRcd" : "cdRc");
    updatePath();
    legend << "\\theta u_x, \\theta u_y, \\theta u_z,";
  }

  simulator.initLog(exp_id, legend3D);

  // plot feature error
  s_.resize(s_rows, false);
  e_.resize(s_rows, false);
  L_.resize(s_rows, 6, false);


  if(e_.size())
  {
    std::string legend_s(legend.str());
    legend_s.back() = ']';
    simulator.logger.save(e_, "err", legend_s, "Feature error");
  }



}
