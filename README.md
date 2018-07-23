# NeoArealize

A model of reaction-diffusion pattern formation in neocortex.

Pre-requisites:

Build and install morphologica from:

https://github.com/ABRG-Models/morphologica

Make sure these packages are installed (Debian/Ubuntu example):

sudo apt install python python-numpy xterm

Now build:

```bash
cd NeoArealize
mkdir build
cd build
cmake ..
make
cd ..
```

To run:

```bash
python sim2.py
```
