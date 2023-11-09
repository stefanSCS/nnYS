# nnYS


Neural networks applied to the modeling of yield surfaces (and functions) in metal plasticity.   

Repository for code and data supporting the article ***"On the use of neural networks in the modeling of yield surfaces"***   

Preprint: [ResearchGate](https://www.researchgate.net/publication/371681181_On_the_use_of_neural_networks_in_the_modeling_of_yield_surfaces)  



There are three directories: **partA**,  **partB** and **partC**.

### **partA** 
Contains all notebooks referenced in Sections 3-6 of the article (direct modeling of yield surfaces via NNs) and also recorded in the Appendix/Table. 


### **partB**   
Contains all notebooks and data supporting Section-7 of the article 
(Mapping the convexity domain of a yield function)

The data files in **partB** are as follows: 

- Training data: ```aCVXbdry_deg_4_test.txt```
- Near boundary points in the training data (See Article/Fig.20 for explanations):
```aCVXbdry_deg_4.txt``` 
- Validation data: ```aCVXbdry_deg_4_Valid_test.txt```  
- Near boundary points in the validation data (See Article/Fig.20 for explanations): ```aCVXbdry_deg_4_Valid.txt``` 


Other files/notebooks in **partB** are as follow:


- NN training:  ```aFitCVXdeg4_12_C_9D4.ipynb```
- NN validation: ```nb24_aValidCVXdeg4_12_C_9D4.ipynb```  
- Poly4 calibration example (used to generate the model in Fig.22): ```nb24_aPolyNoptim.ipynb```
- Supporting python script: ```cvxDomainPolyN.py```  


### **partC**   
(To be released soon) Contains a UMAT implementation of a feed forward (densely connected) NN and some simulation tests. 


### Notes

For convenience (quick preview), I added ```html``` exports of some of the notebooks. Note, however, that the final version of code is always in the corresponding notebook.  


The utility scripts ```ysnnutil.py``` and ```cvxDomainPolyN.py``` are required imports for **partA** and **partB**, respectively. In particular, if Google Colab is used, then variable ```wdir``` must be set to the corresponding path, e.g.    
   ```wdir = '/content/drive/MyDrive/.../partA'```
   
   
   
   