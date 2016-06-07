MODULE PhysicalConstantsModule

  USE KindModule, ONLY: &
    DP

  IMPLICIT NONE
  PRIVATE

  REAL(DP), PUBLIC, PARAMETER :: &
    SpeedOfLightMKS          = 2.99792458e8_DP, &
    GravitationalConstantMKS = 6.673e-11_DP, &
    BoltzmannConstantMKS     = 1.3806503e-23_DP, &
    ElectronVoltMKS          = 1.602176462e-19_DP

END MODULE PhysicalConstantsModule