# Evaluation of efficiency and sensitivity of 1D and 2D sample pooling strategies for SARS-CoV-2 RT-qPCR screening purposes
In this GitHub repository, you can find the code we used for the simulations and the dataset we used for sampling both for the bioRxiv and the manuscript.

## bioRxiv
### _filtered_cq.csv_
This .csv file contains the filtered dataset that was used for the sampling to set the initial situations in the simulations. We applied the following filters to the raw data.
* Full qPCR plates
* Less than 10 positive samples per RNA plate
* Good positive and negative controls

We obtained the correction for the Cq values by calculating the mean Cq of the two positive controls per qPCR plate and calculating the difference of these averages with the global average of the positive controls. 

#### Column information
1. ```plate_id```: Unique identifier of the qPCR plate the sample belongs to.
2. ```sample```: Unique identifier of the sample.
3. ```sample_type```: Since only positive samples are included, the sample_type is always ```std```(for standard).
4. ```target```: The target is always the E gene.
5. ```target_type```: The target type is always the target of interest (```toi```).
6. ```dye```: FAM was used for all qPCR reactions.
7. ```reaction```: Before May 25th, two singleplex reactions were used, after that it was duplex.
8. ```wells```: The well of the reaction. 
9. ```corr_cq```: The Cq value corrected as explained above.

### *sim_code.R*
This is the code we used to perform the simulations. The code is commented and divided in sections. It includes the simulations, data wrangling and code for the figures.


## cid
### _data_first_wave.csv_ and _data_second_wave.csv_
These datasets are only filtered for high-confidence data (good negative and positive controls). They are also not corrected for interplate differences.

### Column information
1. ```sample```: Unique identifier of the sample.
2. ```date```: Date and time at which the corresponding run was initiated.
3. ```cq```: Cq value 

### *simulations.R*
This is the code we used to perform the simulations.

