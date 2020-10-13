#! /bin/bash

echo "|----------------------  gcc --version  ----------------------"
gcc --version

echo "|----------------------  gfortran --version  -----------------"
gfortran --version

echo "|----------------------  mpif90 --version  -------------------"
mpif90 --version

echo "|----------------------  mpif90 -showme  ---------------------"
mpif90 -showme

echo "|----------------------  rm *.o  *.mod *.so *.f90  ------------"

echo "REMOVING ALL *.f90 SOURCE CODES !!!"
rm -I  *.f90

echo "REMOVING ALL *.o , *.mod and *.so FILES !!!"
rm -I  *.o  *.mod  *.so

echo "|----------------------  copying files  ...  -----------------"

# please change folder name
cd ~/abinitio/qe-6.6/Modules/
cp atomic_wfc_mod.f90  output_tau.f90        start_k.f90 \
cell_base.f90         io_global.f90  parallel_include.f90  symm_base.f90 \
cellmd.f90            ions_base.f90  parameters.f90        symmetrize_at.f90 \
constants.f90         kind.f90       parser.f90            util_param.f90 \
control_flags.f90     latgen.f90     random_numbers.f90    volume.f90 \
cryst_to_car.f90      wannier_new.f90  error_handler.f90 \
mp.f90         recips.f90            wypos.f90 \
mp_images.f90  sort.f90 			 invmat.f90 \
~/jianguoyun/Workspace/QE/minimal/src/ > /dev/null 2>&1

# please change folder name
cd ~/abinitio/qe-6.6/UtilXlib/
cp atomic_wfc_mod.f90    output_tau.f90        start_k.f90 \
cell_base.f90         io_global.f90  parallel_include.f90  symm_base.f90 \
cellmd.f90            ions_base.f90  parameters.f90        symmetrize_at.f90 \
constants.f90         kind.f90       parser.f90            util_param.f90 \
control_flags.f90     latgen.f90     random_numbers.f90    volume.f90 \
cryst_to_car.f90      wannier_new.f90  error_handler.f90 \
mp.f90         recips.f90            wypos.f90 \
mp_images.f90  sort.f90 			 invmat.f90 \
~/jianguoyun/Workspace/QE/minimal/src/ > /dev/null 2>&1

# please change folder name
cd ~/abinitio/qe-6.6/PW/src/
cp atomic_wfc_mod.f90    output_tau.f90        start_k.f90 \
cell_base.f90         io_global.f90  parallel_include.f90  symm_base.f90 \
cellmd.f90            ions_base.f90  parameters.f90        symmetrize_at.f90 \
constants.f90         kind.f90       parser.f90            util_param.f90 \
control_flags.f90     latgen.f90     random_numbers.f90    volume.f90 \
cryst_to_car.f90      wannier_new.f90  error_handler.f90 \
mp.f90         recips.f90            wypos.f90 \
mp_images.f90  sort.f90 			 invmat.f90 \
~/jianguoyun/Workspace/QE/minimal/src/ > /dev/null 2>&1

# please change folder name
cd ~/jianguoyun/Workspace/QE/minimal/src/

echo "|----------------------  compile -----------------------------"

GFORT=$(mpif90 -showme)
alias fffc="${GFORT} -std=f2003 -cpp -shared -fPIC -O3 -L/usr/local/lib/ -lblas -llapack -c"
shopt -s expand_aliases

fffc parallel_include.f90
fffc util_param.f90
fffc mp.f90
fffc error_handler.f90

fffc kind.f90
fffc io_global.f90
fffc parameters.f90

fffc mp_images.f90
fffc parser.f90

fffc constants.f90
fffc sort.f90
fffc cryst_to_car.f90
fffc atomic_wfc_mod.f90
fffc control_flags.f90
fffc invmat.f90
fffc latgen.f90
fffc random_numbers.f90
fffc wypos.f90
fffc wannier_new.f90
#fffc input_parameters.f90
#fffc read_cards.f90
fffc cell_base.f90
fffc symm_base.f90
fffc ions_base.f90
fffc output_tau.f90

${GFORT} -std=f2003 -cpp -shared -fPIC -O3  *.o  volume.f90 recips.f90 ./1/main.f90  -L/usr/local/lib/ -lblas -llapack  -o  QE_minimal.so

echo "|----------------------  finalize ----------------------------"
mv QE_minimal.so ../lib/