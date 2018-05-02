The implementations were adapted from their originals in order to receive the
same command line parameters, which are:

Parameter | Type | Description
--- | :---: | ---
```-f``` | Mandatory | The path to the QAP instace file
```-s``` | Optional | The initial seed for the pseudo-random number generator
```-i``` | Optional | The maximum number of iterations to be run
```-o``` | Optional | The output file path

Therefore, an algorithm can be executed with the following command:

	$ ./bin/bls -f ../qap_instances/nug30.dat -s 62912 -i 3000 -o nug30_bls_output.txt
