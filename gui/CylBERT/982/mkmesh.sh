MESH=mesh/mesh
PRIM=${MESH}Prim
mkdir -p mesh
DIA=0.3
HEI=0.6
polyCreateCube -v -Z -s32 -m0 -a0.00001 $MESH
polyTranslate -z -0.5 $MESH
polyScale -x $DIA -y $DIA -z $HEI $MESH
cat data.shm|head -n 82 |tail -n 80 > elec.xyz
polyAddVIP -f elec.xyz -m -99 $MESH
polyAddVIP -x 0 -y 0 -z 0 -m -999 $MESH
polyAddVIP -x 0 -y 0 -z -$HEI -m -1000 $MESH
cp $MESH.poly $PRIM.poly
polyRefineVIPS -r 0.005 $PRIM
tetgen -pazVACq1.2 $PRIM
meshconvert -vDBM -p -it -o $PRIM $PRIM.1
rm -f $PRIM.1.*
tetgen -pazVACq1.2 $MESH
meshconvert -vDBMV -it -o $MESH $MESH.1
rm -f $MESH.1.*

