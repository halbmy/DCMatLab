sh mkmesh.sh
dcmod -v -1 -s data.shm -R -m mesh/mesh.bms -p primpot mesh/meshPrim.bms
dcedit -v -D -c data.ohm -f "a b m n k" -o data.def data.ohm
