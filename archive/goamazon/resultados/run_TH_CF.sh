set -e
cd /mnt/c/git_repos/SingleColumnModelMicrophysicsAndLightning
test "$(realpath build)" = /mnt/c/git_repos/SingleColumnModelMicrophysicsAndLightning/build
backup=$(mktemp -d resultados/pre_TH_CF_XXXXXX)
cp build/exec_TH_CF.log build/namelist.mp_offline "$backup/"
make clean > resultados/build_TH_CF.log 2>&1
make -j4 >> resultados/build_TH_CF.log 2>&1
cp namelist.mp_offline build/namelist.mp_offline
sed -i 's|data/GOAMAZON_SINGLE_SCM_driver.nc|../data/GOAMAZON_SINGLE_SCM_driver.nc|' build/namelist.mp_offline
cd build
./mp_offline_driver > exec_TH_CF.log 2>&1
ls -lh output_mp_offline.nc
tail -n 25 exec_TH_CF.log
ncdump -h output_mp_offline.nc > ../resultados/TH_CF_header.txt
