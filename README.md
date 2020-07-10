# Evaluation of efficiency and practicality of 1D and 2D sample pooling strategies for diagnostic screening purposes
In this GitHub repository, you can find the code we used for the simulations and the dataset we used for sampling.

## filtered_cq.csv
This .csv file contains the filtered dataset that was used for the sampling to set the initial situations in the simulations. We applied the following filters to the raw data.
*Full qPCR plates
*Less than 10 positive samples per RNA plate
*Good positive and negative controls

We obtained the correction for the Cq values by calculating the mean Cq of the two positive controls per qPCR plate and calculating the difference of these averages with the global average of the positive controls. 
