# ImagCNV
a novel framework integrated SNV-based recalibration probabilistic model and image classification architecture for CNVs discovery

`Because of the large volume of training datasets, click the following link to get the whole scripts and data.`  
[Google Drive](https://drive.google.com/file/d/16-QvI2g-TIMFEPOf3VCjGRr3Hujl-Rcz/view?usp=sharing)

## Features
* Integration of RD, PE+SR approaches, and optimize bin size of RD approach
* Quality-control strategies
* A Naive Bayesian Model to recalculate the CNV probability 
* Image classification by deep neural network

## Prerequisites
* Linux or UNIX
* Python >=2.7.10 (For Quality-control strategies and Naive Bayesian Model)
* Python3 >=3.6.3 (For TensorFlow only)
* [Tensorflow](https://tensorflow.google.cn/install) (pip3 install tensorflow)
* CNVnator >=v0.2.7
* LUMPY >=v 0.2.7
 
## Processing flow
### 1. Run the Quality-control strategies for RD and PE+SR methods
```
python bin/RD_PE_SR_optimal.py test_sample
```
The input files included the output of CNVnators and LUMPY, `test_sample.cnv.xls` and `test_sample.sv.result.xls` under the directory of `example`.  
It will generate `test_sample.merged.final.test.xls` after optimal Quality-control strategies.
 
### 2. Image Classification
```
python  bin/InceptionV3/image_identify.py  example/test_sample.image_classify_input.xls  example/pic  bin/InceptionV3/  example/  test_sample
```
* Input file: `example/test_sample.image_classify_input.xls`  
* the images of CNV depth corresponding to input file: `example/pic`  
* the retrained the last layer of InceptionV3: `bin/InceptionV3/`  
* output directorty: `example/`  
* sample name: `test_sample`  

Then generate the output file `test_sample.image.fil.xls` under the `example` directory
  
You can retrain the last layer of InceptionV3 ownself using the code below:
```
python3 bin/InceptionV3/retrain.py --image_dir bin/InceptionV3/pic_training/ --bottleneck_dir bin/InceptionV3/bottleneck/ --how_many_training_steps 200 --model_dir bin/InceptionV3/model_classify/ --output_graph bin/InceptionV3/output_graph.pb --output_labels bin/InceptionV3/output_labels.txt --saved_model_dir bin/InceptionV3/model_out
  ```
  
### 3. Implement the Naive Bayesian model
```
python bin/naive_bayes.py  example/bayes/CNV_input/  example/bayes/training_datasets/normal/  example/bayes/training_datasets/delete/  example/bayes/variants/  example/bayes/output.bayes.cnv.xls 50
```
* input files of CNVs: `example/bayes/CNV_input/`  
* the training datasets of normal regions: `example/bayes/training_datasets/normal/`  
* the training datasets of CNV regions: `example/bayes/training_datasets/delete/`  
* the SNP/InDel variants: `example/bayes/variants/`  
* output file: ` example/bayes/output.bayes.cnv.xls`  
* the shielding parameter: `50`

The output format of the Comprehensive framework as follows:
```
chr1    1011001     1014000    deletion        0       1
chr1    1285001     1287000    deletion        0       1
chr1    17185001    17277000   duplication     0       0
chr1    64839001    64852000   deletion        1       1
```
The first four columns are the putative regions of CNVs calling by RD,PE+SR approaches. The fifth column is the predicition of Image Classification, and `0` represents NO and `1` YES. The last columns is the result of naive bayesian model. The candidate CNVs are marked `1` in the fifth column or the last column.


## Performace
The novel framework integrated SNV-based recalibration probabilistic model and image classification architecture to infer CNVs can achieve significantly lower FDR and comparative sensitivity than the separate RD model, or PE and SR models. The multiple modules used in the framework can emphasis on different key points which increase the robustness of the framework.



