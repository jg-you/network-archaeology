# Compile OD
cd tools/bins/
echo "[Entering $PWD]"
echo "[Compiling OD.cpp]"
g++ -o3 OD.cpp -o OD
echo "[Compilation of OD.cpp done]"
# Go back to root
cd ../../
echo "[Entering $PWD]"
# Compile generator
cd generators
echo "[Entering $PWD]"
echo "[Launching compilation of the network generators (CMAKE)]"
cmake .;
make
echo "[Compilation of the network generators done]"
cd ..
echo "[Entering $PWD]"
echo "[Done]"