# HeadTailAssigner (HTA) Non-Interactive version

Welcome to the HTA non-interactive version, where you can learn how to use HTA to run head and tail assignment in batch in a terminal.

## Documentation

### How to run HTA

1. After downloading the package in a desired folder, create the conda environment as explained in the main [README.md](README.md)

2. Use the [input.csv](examples\non-interactive\input.csv) provided in this folder to run HTA. If you would like to prepare your own input, please follow the same format as the example input.

3. In the terminal, move to the folder where main.py is located. Be aware that you always have to be in the same folder as the main.py file, since python needs to find this file to run HTA.

4. To run the code, type at your terminal:

```
python main.py examples\non-interactive\input.csv
```

5. Now HTA is going to run and you can see at your terminal the entire process being logged. If you would like to check this information later, this log is saved in the output-directory folder as log.out.

6. When the run finishes, you can check output-directory folder to verify all the files that were created. HTA creates an image of the monomers with and without index in png, the input for GAMESS in a inp file, the log of the GAMESS run in log file, the structure of the monomer in xyz file and also the files with the extracted Mulliken's population information in txt and the indexes that present the highest electronic population in csv file. So you can check how HTA produces the head and tail assignment using these files.

7. The results of the head and tail assignment can be located at the results.csv file in the same directory as the main.py file.

---

## Observations

- If you have any questions be welcomed to contact us by opening an issue on the [HTA github](https://github.com/IBM/HeadTailAssign). We will try to answer as soon as possible. If you use HTA in your research please do not forget to cite us.

---

## Authorship

- Author: **Brenda Ferrari** [brendaferrari](https://github.com/brendaferrari)
- Co-author: **Ronaldo Giro** [rgiro](https://github.com/rgiro)
