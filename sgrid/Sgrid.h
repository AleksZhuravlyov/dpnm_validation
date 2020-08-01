#include<string>
#include<vector>
#include<map>

#include <Eigen/Dense>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>


class Sgrid {

public:

  Sgrid(Eigen::Ref<Eigen::Vector3i> pointsDims,
        Eigen::Ref<Eigen::Vector3d> pointsOrigin,
        Eigen::Ref<Eigen::Vector3d> spacing);

  Sgrid(const std::string &fileName);

  virtual ~Sgrid() = default;

  /// Main methods

  void save(const std::string &fileName);


  /// Accessors and mutators

  Eigen::Ref<Eigen::Vector3i> getPointsDims();

  void setPointsDims(Eigen::Ref<Eigen::Vector3i> pointsDims);


  Eigen::Ref<Eigen::Vector3d> getPointsOrigin();

  void setPointsOrigin(Eigen::Ref<Eigen::Vector3d> pointsOrigin);


  Eigen::Ref<Eigen::Vector3d> getSpacing();

  void setSpacing(Eigen::Ref<Eigen::Vector3d> spacing);


  Eigen::Ref<Eigen::Vector3i> getCellsDims();

  void setCellsDims(Eigen::Ref<Eigen::Vector3i> cellsDims);


  std::map<int, Eigen::Ref<Eigen::Vector3i>> getFacesDims();

  void setFacesDims(std::map<int, Eigen::Ref<Eigen::Vector3i>> facesDims);


  Eigen::Ref<Eigen::Vector3d> getFaceS();

  void setFaceS(Eigen::Ref<Eigen::Vector3d> faceS);


  std::map<int, Eigen::Ref<Eigen::VectorXi>> getNeighborsFaces();

  void setNeighborsFaces(std::map<int, Eigen::Ref<Eigen::VectorXi>> neighborsFaces);

  std::map<int, Eigen::Ref<Eigen::VectorXi>> getNeighborsCells();

  void setNeighborsCells(std::map<int, Eigen::Ref<Eigen::VectorXi>> neighborsCells);


  std::map<int, Eigen::Ref<Eigen::VectorXi>> getNormalsNeighborsCells();

  void setNormalsNeighborsCells(
      std::map<int, Eigen::Ref<Eigen::VectorXi>> normalsNeighborsCells);

  std::map<int, Eigen::Ref<Eigen::VectorXi>> getNormalsNeighborsFaces();

  void setNormalsNeighborsFaces(
      std::map<int, Eigen::Ref<Eigen::VectorXi>> normalsNeighborsFaces);


  std::map<std::string, Eigen::Ref<Eigen::VectorXd>> getPointsArrays();

  void setPointsArrays(
      std::map<std::string, Eigen::Ref<Eigen::VectorXd>> pointsArrays);

  std::map<std::string, Eigen::Ref<Eigen::VectorXd>> getCellsArrays();

  void setCellsArrays(
      std::map<std::string, Eigen::Ref<Eigen::VectorXd>> cellsArrays);

  std::map<std::string, Eigen::Ref<Eigen::VectorXd>> getFacesArrays();

  void setFacesArrays(
      std::map<std::string, Eigen::Ref<Eigen::VectorXd>> facesArrays);


  Eigen::Map<Eigen::Vector3i> _pointsDims;
  Eigen::Map<Eigen::Vector3d> _spacing;
  Eigen::Map<Eigen::Vector3d> _pointsOrigin;

  int _pointsN;

  Eigen::Map<Eigen::Vector3i> _cellsDims;
  int _cellsN;

  std::map<int, Eigen::Map<Eigen::Vector3i>> _facesDims;
  int _facesN;

  double _cellV;
  Eigen::Vector3d _faceS;

  std::map<int, Eigen::Map<Eigen::VectorXi>> _neighborsFaces;
  std::map<int, Eigen::Map<Eigen::VectorXi>> _neighborsCells;

  std::map<int, Eigen::Map<Eigen::VectorXi>> _normalsNeighborsCells;
  std::map<int, Eigen::Map<Eigen::VectorXi>> _normalsNeighborsFaces;


  std::map<std::string, Eigen::Map<Eigen::VectorXd>> _pointsArrays;
  std::map<std::string, Eigen::Map<Eigen::VectorXd>> _cellsArrays;
  std::map<std::string, Eigen::Map<Eigen::VectorXd>> _facesArrays;


  /// Accessory shitty constructor methods

  void calculatePointsN();

  void calculateCellsDims();

  void calculateCellsN();

  void calculateFacesDims();

  void calculateFacesN();

  void calculateCellV();

  void calculateFaceS();

  void calculateNeighborsFaces();

  void calculateNeighborsCells();

  void calculateNormalsNeighborsCells();

  void calculateNormalsNeighborsFaces();


  void calculateGridProps();


};


void saveFilesCollectionToFile(const std::string &fileName,
                               const std::vector<std::string> &filesNames,
                               const std::vector<std::string> &filesDescriptions);