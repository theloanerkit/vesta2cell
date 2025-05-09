# vesta2cell
Converts a [VESTA](https://jp-minerals.org/vesta/en/) ```.vesta``` output file into a [CASTEP](https://www.castep.org/) ```.cell``` file.

### Using the command line
To run using the command line:
```
python3 vesta2cell.py -s [SEED] -o [OUTPUT] -sp [nospin/collinear/noncollinear]
```
- ```-s``` specifies the seed of the input ```.vesta``` file (required)
- ```-o``` specifies the seed of the output ```.cell``` file (optional)
    - if not specified, the input seed will be used for the ```.cell``` file
- ```-sp``` specifies the treatment of spin when converting (optional)
    - default is ```nospin```

### From within another python file
```
import vesta2cell as vc

vc.convert([SEED],[OUPTUT],[nospin/collinear/noncollinear])
```