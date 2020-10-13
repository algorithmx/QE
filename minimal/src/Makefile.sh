#! /bin/bash

QE_source_root_fd="/home/dabajabaza/abinitio/qe-6.6"
Workspace="/home/dabajabaza/jianguoyun/Workspace/QE/minimal/src"

fd_Modules=(\
"kind.f90" \
"constants.f90" \
"parameters.f90" \
"io_global.f90" \
"cell_base.f90" \
"ions_base.f90" \
"recips.f90" \
"volume.f90" \
"invmat.f90" \
"wypos.f90" \
"wannier_new.f90" \
"control_flags.f90" \
"latgen.f90" \
"sort.f90" \
"random_numbers.f90" \
"parser.f90" \
"mp_images.f90" \
)

fd_PW_src=(\
"atomic_wfc_mod.f90" \
"symm_base.f90" \
"symmetrize_at.f90" \
"output_tau.f90" \
)

fd_UtilXlib=(\
"parallel_include.f90" \
"util_param.f90" \
"error_handler.f90" \
"mp.f90" \
)




echo "|----------------------  gcc --version  ----------------------"
gcc --version

echo "|----------------------  gfortran --version  -----------------"
gfortran --version

echo "|----------------------  mpif90 --version  -------------------"
mpif90 --version

echo "|----------------------  mpif90 -showme  ---------------------"
mpif90 -showme

echo "|----------------------  rm *.o  *.mod *.so *.f90  ------------"

echo "| REMOVING ALL *.f90 SOURCE CODES !!!"
rm -I  *.f90

echo "| REMOVING ALL *.o , *.mod and *.so FILES !!!"
rm -I  *.o  *.mod  *.so

echo "|----------------------  copying files  ...  -----------------"

for f in "${fd_Modules[@]}" 
do
    cp $QE_source_root_fd/Modules/$f $Workspace > /dev/null 2>&1
done

for f in "${fd_UtilXlib[@]}" 
do
    cp $QE_source_root_fd/UtilXlib/$f $Workspace > /dev/null 2>&1
done

for f in "${fd_PW_src[@]}" 
do
    cp $QE_source_root_fd/PW/src/$f $Workspace > /dev/null 2>&1
done

echo "| Done."

echo "| ----------------------  compile -----------------------------"

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
#fffc cryst_to_car.f90
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

echo "| Done."

echo "| ----------------------  finalize ----------------------------"

mv QE_minimal.so ../lib/

echo "| Done."