# Instructions for settings and analysis
ALC_TRAJECTORY requires two files:

* TRAJECTORY: must contain the atomic positions recorded along the MD trajectory. To date, formats "xyz" and "vasp" (XDATCAR) are the only two implemented options.
* SETTINGS: contains the instructions for data analysis.

The SETTINGS file must be present in the folder where the ALC_TRAJECTORY is executed from, otherwise the program will print an error message and abort. Comments to this file can be added using the symbol "#". It is recommended to add a descriptive header in the first lines of the file for revision purposes. Different files will be printed depending on the selected options. The code is flexible enough to recognize directives independently of capitalization. For example, directive ***Ensemble***, ***eNsemBle***, ***ENSEMBLE***, etc, will all be interpreted as ***ensemble***. Together with the individual directives, it is also required to specify blocks, which are declared with the character "***&***", followed by the ***name*** of the block (i.e. ***&name***), and must be closed with ***&end_name***. ALC_TRAJECTORY generates the OUTPUT file, which details the input settings (based on the content of the SETTINGS file) and the generated information from the data analysis.

Upon execution, the code first checks the correctness of the syntax and format for the defined settings. If a problem is found, an error message is printed to OUTPUT file and the execution is aborted. The code then reads the TRAJECTORY file and the collected information is compared against the directives of the SETTINGS file. In case there is an inconsistency but the calculation can still proceed, ALC_TRAJECTORY prints a warning message. If the specified settings are incorrect, the program will print an error message instructing the user what to fix (hopefully). 

The structure of the SETTINGS file can be divided in model and trajectory related directives/blocks. In the following sections, we will use the example case of nano-confined water in Nafion to explain the different functionalities. This model system constitutes an example of a reactive system, where chemical species change along the trajectory. The SETTINGS file can be found in the ***examples/Nafion-hydrated/*** folder. The TRAJECTORY file is also attached, so the user can run the code and change with the settings at will.  

The files for the analysis of bulk liquid water (modelled by 64 water molecules in a cubic box) can be found in the ***examples/bulk-water/*** folder. In contrast to the Nafion, the model for bulk liquid water constitutes an example of a **non-reactive** system. We shall not discuss this case in the following notes.  In contrast to reactive systems, the definition of the ***&search_chemistry*** block  (see the section "Identification of chemistry changes" below) is NOT NEEDED for non-reactive systems.   

Accepted units are fs (femtoseconds) and ps (picoseconds) for time-related directives, and Angstrom and Bohr for distance-related directives (cutoff). Before proceeding to the desciption of settings and funcionalities, we define the following acronyms:  

* **TCF**: Transfer Correlation Function  
* **OCF**: Orientational Correlation Function  
* **RDF**: Radial Distribution Function  
* **MSD**: Mean Square Displacement  

## <span style="color: green">  Model settings</span>
To describe the implemented functionalities for analysis, we will use the example case of nano-confined water in Nafion. The SETTINGS file can be found in the ***examples/Nafion-hydrated/*** folder.

### <span style="color: maroon">Format for the atomic positions</span>
The format of the the TRAJECTORY file will be determined by the computational code used to run the MD simulation. To date, ALC_TRAJECTORY only accepts "xyz" and "vasp" formats. The format for the geometry of the trajectory is specified as follows:  

***input_geometry_format***&nbsp;&nbsp;xyz

### <span style="color: maroon"> Input composition</span>
To specify the atomic details of the input model, it is necessary to define the ***&input_composition*** block. This block allows tagging the atoms of the model according to the structure of the system under consideration. Tags are needed to classify elements with different chemical environments. For example, a model can have Hydrogen atoms as part of water molecules and Hydrogen atoms as part of a modelled backbone for a membrane: the chemical element is the same (H), but the chemical environment around each hydrogen is completely different. The structure of this block for the present example is defined as follows:  

***&input_composition***  
&nbsp;&nbsp;&nbsp;atomic_species&nbsp;&nbsp;  7  
&nbsp;&nbsp;&nbsp;tags     &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Cch&nbsp; Hch&nbsp;  Hs&nbsp;  Och&nbsp;  Sch&nbsp; Hw&nbsp;&nbsp; Ow  
&nbsp;&nbsp;&nbsp;amounts  &nbsp; 56&nbsp;&nbsp;&nbsp;       60&nbsp;&nbsp;&nbsp;  4&nbsp;&nbsp;&nbsp;&nbsp;  12 &nbsp;&nbsp;           4 &nbsp;&nbsp;       112 &nbsp; 56  
&nbsp;&nbsp;      elements &nbsp;  C&nbsp;&nbsp;&nbsp;&nbsp;  H &nbsp;&nbsp;&nbsp;       H&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;  O&nbsp;&nbsp;&nbsp;&nbsp; S &nbsp;&nbsp;&nbsp;   H &nbsp;&nbsp;&nbsp;&nbsp; O  
***&end_input_composition***

The first directive must be ***atomic_species***, in this case equal to 7. This setting indicates that 7 different atomic tags will be used to label the atoms. The next task is to define the ***tags***, ***amounts*** and ***elements*** directives. The order for the definition of such directives is irrelevant. These settings describe what type of atoms and how many of them constitute each frame of the trajectory. Atomic tags are selected based on a reference configuration (see the [DL_FIELD](https://pubs.acs.org/doi/10.1021/acs.jcim.6b00323) code for example). In reactive systems, the tags of certain atoms are expected to change with respect to the tagging of the reference configuration. We shall return to this point later when we describe the ***&search_chemistry*** block. Consistent with the definition of the ***&input_composition*** block, each frame (atomic configuration) of the trajectory must have the following structure/order:  

* 56&nbsp;&nbsp; C atoms with tag Cch, followed by  
* 60&nbsp;&nbsp; H atoms with tag Hch, followed by  
* 4&nbsp;&nbsp;&nbsp;&nbsp;  H atoms with tag Hs,  followed by  
* 12&nbsp;&nbsp; O atoms with tag Och, followed by  
* 4&nbsp;&nbsp;&nbsp;&nbsp; S atoms with tag Sch, followed by 
* 112 H atoms with tag Hw, followed by
* 56&nbsp;&nbsp; O atoms with  tag Ow.

The suffixes "ch", "s" and "w" here indicate "chain", "sulphonic" and "water", respectively, but this choice is arbitrary and subject to user's criteria. Ideally, the initial configuration of the MD run should arrange atoms in groups to facilitate the definition of the block. Unfortunately this is not always possible. Still, the block is sufficiently flexible to allow multiple definitions of the same tags. For example, let's assume that the atomic configuration of 30 Hch atoms (out of 60) have been grouped at the end. The structure for the block would be:
 
***&input_composition***  
&nbsp;&nbsp;&nbsp;atomic_species&nbsp;&nbsp;  <span style="color: blue">8</span>  
&nbsp;&nbsp;&nbsp;tags     &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Cch&nbsp; Hch&nbsp;  Hs&nbsp;  Och&nbsp;  Sch&nbsp; Hw&nbsp;&nbsp; Ow &nbsp;&nbsp; <span style="color: blue">Hch</span>  
&nbsp;&nbsp;&nbsp;amounts  &nbsp; 56&nbsp;&nbsp;&nbsp; <span style="color: blue">30</span>&nbsp;&nbsp;&nbsp;  4&nbsp;&nbsp;&nbsp;&nbsp;  12 &nbsp;&nbsp;           4 &nbsp;&nbsp;       112 &nbsp; 56&nbsp;&nbsp; <span style="color: blue">30</span>  
&nbsp;&nbsp;      elements &nbsp;  C&nbsp;&nbsp;&nbsp;&nbsp;  H &nbsp;&nbsp;&nbsp;       H&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;  O&nbsp;&nbsp;&nbsp;&nbsp; S &nbsp;&nbsp;&nbsp;   H &nbsp;&nbsp;&nbsp;&nbsp; O&nbsp;&nbsp;&nbsp; <span style="color: blue">H</span>  
***&end_input_composition***

It is important to clarify that although the last group of 30 Hch atoms are the same as the previous 30 Hch, it is needed to increase "atomic_species" by 1, from 7 to 8. If there is an inconsistency between the settings in the block and the structure of the TRAJECTORY file, the code will inform the user and abort the execution.

### <span style="color: maroon"> Simulation cell</span>
Only for those trajectories recorded in "xyz" format (see the ***input_geometry_format*** directive above), the code will ask to define the simulation cell, which remains unchanged along the trajectory. In fact, trajectories in "xyz" format are only compatible with NVE and NVT simulations (see the ***ensemble*** directive below). Information for the simulation cell for the current example is set as follows:  

***cell_units*** &nbsp;&nbsp;  Angstrom  

***&simulation_cell***  
&nbsp;&nbsp;             10.064000&nbsp;&nbsp;  0.0000000&nbsp;&nbsp;    0.0000000  
&nbsp;&nbsp;&nbsp;&nbsp;  0.000000             13.0735195&nbsp;&nbsp;    0.0000000  
&nbsp;&nbsp;&nbsp;&nbsp;  0.000000&nbsp;&nbsp;  0.0000000               20.0000000  
***&end_simulation_cell***

Both ***cell_units*** and ***&simulation_cell*** are NOT NEEDED if the TRAJECTORY was recorded in "vasp" format. If there is an inconsistency between this block and the atomic positions of the trajectory, the code will inform the user that something incorrect and abort the execution.

### <span style="color: maroon"> Identification of chemistry changes</span>
Reactive systems are characterized by the breaking and formation of chemical bonds along the trajectory, which leads to the changes of the constituents species. To compute reactive systems, it is first necessary to set 

***change_chemistry***&nbsp;&nbsp;&nbsp;True

which instructs ALC_TRAJECTORY to search for changes of chemical species within the model. Changes are monitored through the identification of donor and acceptor sites. The user needs to define following block: 

***&search_chemistry***  
&nbsp;&nbsp;&nbsp;total_number &nbsp;&nbsp;4  

&nbsp;&nbsp;&nbsp;&bonding_criteria  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;number_of_bonds&nbsp;&nbsp;   3  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;only_element   &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;   H  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;cutoff&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;1.3  &nbsp;&nbsp;Angstrom  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span style="color: blue"> ***&possible_extra_bonds***  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;types_of_bonds&nbsp;&nbsp; 2  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Och&nbsp;&nbsp;Hw&nbsp;&nbsp;1.25&nbsp;&nbsp;  Angstrom  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Och&nbsp;&nbsp;Hs&nbsp;&nbsp;&nbsp;1.25&nbsp;&nbsp;&nbsp;Angstrom  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;***&end_possible_extra_bonds***</span>  
&nbsp;&nbsp;&nbsp;&end_bonding_criteria  

&nbsp;&nbsp;&nbsp;&acceptor_criteria  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;include_tags&nbsp;&nbsp;&nbsp; 2&nbsp;&nbsp;Ow&nbsp; Och  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;exclude_pairs&nbsp;            1&nbsp; Och  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;cutoff       &nbsp;&nbsp; 3.20 &nbsp;&nbsp; Angstrom  
&nbsp;&nbsp;&nbsp;&end_acceptor_criteria  

***&end_search_chemistry***

The ***total_number*** directive sets to identify a total of 4 donor sites (thence chemical species) at each frame of the trajectory. In this example, the total number of 4 comes from the H atom of the 4 SO<sub>3</sub>H groups. Donors constitute reference atomic sites, part of the formed chemical species. For the present case, we set to identify H<sub>3</sub>O<sup>+</sup> species using the ***&bonding_criteria*** sub-block, which specifies that each reference sites (O atoms in this case, see ***&acceptor_criteria***) forms 3 bonds only with H atoms, independently of theirs tag. A bond between the reference site (donor) and the H atoms is subject to a cutoff distance criterion of 1.3 Angstrom, beyond which the bond is considered to be broken.  
This information would be enough to track H<sub>3</sub>O<sup>+</sup> in liquid phase. For this example of Nafion, however, H<sub>2</sub>O and H<sub>3</sub>O<sup>+</sup> species are not part of bulk liquid water but constitute an interface with the backbone structure of the membrane. Consequently, protons can also form bonds with the oxygen atoms of the SO<sub>3</sub></sup>-</sup> groups. Such oxygens are labelled as Och in the ***input_composition*** block. To account for these bonds, the user must define the ***&possible_extra_bonds*** sub-block and specify the number of type of bonds (2 in this case). Together with the atomic tags the define the bond, the user must also define the cutoff value and the units as shown. This ends the specifications for the ***&bonding_criteria*** block. To indicate that the possible tags (Ow and Och) become the reference site of a newly formed chemical species, the implemented algorithm retags the sites by adding an asterisk symbol `*` as apex, so the tags become Ow`*` and Och`*`. The same is done with the corresponding H atoms. This retagging is important when defining the settings for RDF analysis.  

So far, we have implicitly discussed about the donors (Ow and Och) based on chemical knowledge for the system under consideration. This is not enough for executing the code. In fact, we need to specify the criteria not only to identify the donor sites but also to track the transferring of atomic species to neighbouring acceptors. This information must be added using the ***&acceptor_criteria*** block. The ***include_tags*** directive sets the number of atomic tags (2, Ow and Och) that are possible candidates to become the reference site of the chemical species. The implemented algorithm first identifies the reference species sites, which become the donors. For the subsequent frames of the MD trajectory, ALC_TRAJECTORY determines the bonding pattern for all those Ow and Och neighbouring sites, only within the environment region around each donors. The environment is defined using a distance cutoff criterion, the ***cutoff*** directive, of 3.2 Angstrom (note: this cutoff is different from the bonding cutoff. A reasonable goof value for this directive can be obtained from the RDF analysis). Evaluating the bonding criteria for donors and possible acceptors allows detecting (proton) transfer events and tracking changes of the chemical species. The settings for the ***include_tags*** directive also indicate that when the site is Ow`*` or Och`*`, both Ow and Och can be considered as acceptors within the distance cutoff. For this system, we know that Och`*`-> Och transitions are unphysical. To accelerate the search, we can (optionally) opt to discard the calculation of Och-Och pairs in the search. This is specified with the ***exclude_pairs*** as shown.  
The ***&search_chemistry*** block prints the TRACK_CHEMISTRY file with the xyz position for the donor sites along the trajectory. As explained, such donors can be interpreted as the location of the chemical species, in this case H<sub>3</sub>O<sup>+</sup> and SO<sub>3</sub>H. The residence times are reported in the RES_TIMES file for each of the species, together with the corresponding tag (either Ow`*` or Och`*`). The percentage of residence for the tags, accounting for all the 4 species, are printed as a table in the OUTPUT file. Finally, if the ***&lifetime*** block is defined (see below), ALC_TRAJECTORY also compute the Transfer Correlation Function (TCP), being the proton-TCF in the present case.  
It is important to emphasise that the setting of the ***&search_chemistry*** block ONLY requires atomic labels, defined in ***&input_composition***, and the definition of few cutoff values. This structure for the settings is rather flexible and general. Most importantly, it avoids the complexity of setting system-specific bond lengths and angle parameters, as required by other analysis tools.

### <span style="color: maroon"> Definition of monitored species</span>
The model related settings of the previous sections are already sufficient to evaluate chemistry changes along the trajectory.  However, it is also important to evaluate species that are not the chemical active species, but part of the environment (solution). This includes the possibility to perform Mean Square Displacement (MSD) and  Orientational Correlation Function (OCF) analyses for the non-reactive species of the system. We shall define and refer to such species as "monitored species". In the present case, monitored species are the nano-confined water molecules. The definition of the monitored species for this case is trivial and straightforward, but illustrates the flexibility to define (other) species:  

***&monitored_species***  
&nbsp;&nbsp;&nbsp;name H2O  
&nbsp;&nbsp;&nbsp;reference_tag Ow  
&nbsp;&nbsp;&nbsp;<span style="color: blue"> ***&atomic_components***  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;number_components 2  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;H &nbsp;&nbsp; 2  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;O &nbsp;&nbsp; 1  
&nbsp;&nbsp;&nbsp;***&end_atomic_components***</span>  
&nbsp;&nbsp;&nbsp;bond_cutoff&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;  1.2 &nbsp;&nbsp;Angstrom  
&nbsp;&nbsp;&nbsp;compute_amount&nbsp;&nbsp;  .True.  
***&end_monitored_species***

In this example, we name the species as" H2O" and use the Ow as the reference atomic tag. The atomic composition is set using the ***&atomic_components*** block, obvious in this case. The cutoff sets that the maximum bonding distance criterion between the O and the H atoms for the set to be considered as a "monitored species". When setting the optional ***compute_amount***  directive to .True. (set to .False. by default), the code computes the average number (and the standard deviation) of monitored species, and the result is printed to the OUTPUT file. The computed number depends on the settings for the &region block and the ignore_initial directive, if defined.

## <span style="color: green"> Trajectory settings for analysis</span>
Having defined the model related settings, the next step is to specify trajectory related directives and blocks for MD analysis. The first two compulsory directives must be ensemble used and the timestep to record the trajectory. It is important to clarify that this timestep is not the time step used for the numerical integration of the equation of motion, but the time step used to record the configurations of the trajectory. For this example we have:  

***ensemble***&nbsp;&nbsp;&nbsp;       NVE  
***timestep***&nbsp;&nbsp;&nbsp;&nbsp; 10&nbsp;&nbsp;  fs  

These two settings are important to evaluate and process the data. It is user's responsibility to ensure the correctness of these directives and avoid an incorrect analysis. Options for ***ensemble*** are NVE, NVT and NPT. Analysis for NPT ensembles is NOT compatible for "xyz" formats because the simulation cell is NOT fixed along the trajectory. In contrast, NPT is a valid option for trajectories in "vasp" format.  
Optionally, for non-reactive systems, the user can request printing a trajectory that contains the labels used to identify the species:
  
***print_retagged_trajectory***&nbsp;&nbsp;&nbsp;TRUE &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; (only if ***change_chemistry***&nbsp;&nbsp;&nbsp;True)

This last setting prints the relabelled trajectory to the TAGGED_TRAJECTORY file, which is useful for debugging purposes.

### <span style="color: maroon">Radial Distribution Function (RDF)</span>
This capability is a common feature in the majority of available software for MD data analysis. The user should use the tags defined in the ***&input_composition*** block to define species of **type a** and **type b**. Tags for both **a** and **b** species should correspond to the same chemical element. The chemical element for type **a** can be different from the chemical elements for the type **b** . In addition, the user needs to specify the discretization for the distance, ***dr***, set to 0.02 Angstrom in this case.

***&rdf***  
&nbsp;&nbsp;&nbsp;dr &nbsp;&nbsp;0.02&nbsp;&nbsp;  Angstrom  
&nbsp;&nbsp;&nbsp;tags_species_a&nbsp;&nbsp;  1&nbsp;&nbsp;  Ow  
&nbsp;&nbsp;&nbsp;tags_species_b&nbsp;&nbsp;  2&nbsp;&nbsp;  Ow* &nbsp;&nbsp; Och*  
***&end_rdf***  

In this block, ALC_TRAJECTORY is instructed to compute the RDF between the water oxygens (Ow) and all the oxygen sites of the chemically formed species (Ow* and Och*). Note the values of 1 and 2, preceding the definition of the tags for type a and b, define the amount of atomic tags to be considered. Results are printed to the RDF file.

### <span style="color: maroon"> Definition of time intervals for data analysis</span>
The computation of TCFs, OCFs and MSDs requires dividing of the whole trajectory in time segments. In ALC_TRAJECTORY, this is done through the definition of the following block, which in this case reads:

***&data_analysis***  
&nbsp;&nbsp;&nbsp;time_interval&nbsp;&nbsp;   10.0&nbsp;&nbsp;  ps  
&nbsp;&nbsp;&nbsp;overlap_time &nbsp;&nbsp;   0.50&nbsp; ps  
&nbsp;&nbsp;&nbsp;ignore_initial&nbsp;&nbsp;&nbsp;  2.00&nbsp;&nbsp;ps  
***&end_data_analysis***  

The ***time_interval*** directive specifies the duration of the time segment used to compute the properties. In this case, the time segment is set to 10 ps. It might well happen that the computation of TCFs, OCFs and/or MSDs require different time intervals, for which the user should compute each quantity separately in different runs.  
Optionally, the user can delay the start of the analysis by using the ***ignore_initial*** directive, which ignores the first part of the trajectory (2 ps here).  
Generally, the size for systems compatible with standard DFT simulations are often not large enough statistically-wise. In addition, the length of computed trajectories are limited to tens of picoseconds. Consequently, the computed properties are subject to significantly large errors. To reduce this effect, we follow the strategy of S. Kim et al. [[1]](https://chemistry-europe.onlinelibrary.wiley.com/doi/10.1002/cphc.202000498) and use multiple time origins, which are separated by the settings of the ***overlap_time*** directive (0.50 ps). In this example the analysis for the first segment starts at 2 ps and it lasts for 10 ps (up to the first 12 ps of the trajectory). At 2.50 ps, a new analysis starts, which also lasts for 10 ps and extends to the 12.50 ps of the trajectory.  Likewise, a third analysis starts at 3.0 ps and finishes at 13 ps. The process is repeated until the starting time of the last segment, which is the total time recorded for the trajectory minus 10 ps. At the end of the cycle, we will have multiple analysis, each of 10 ps length. ALC_TRAJECTORY will also compute the average quantity (AVG), together with the standard deviation (STD), the (AVG+STD) and the (AVG-STD) along the time interval. Results are printed to files TCF_AVG, OCF_AVG and MSD_AVG, depending on which block are defined (see below).

### <span style="color: maroon">Transfer Correlation Function (TCF) </span>
TCFs can be computed only for reactive systems. Two directives are required for this block:

***&lifetime***  
&nbsp;&nbsp;&nbsp;method &nbsp;&nbsp; HiCF  
&nbsp;&nbsp;&nbsp;rattling_wait&nbsp;&nbsp;  0.2&nbsp;&nbsp; ps  
***&end_lifetime***  

The ***method*** directive selects the approximation to compute the TCF. In ALC_TRAJECTORY, both History-Independent and History-Dependent Correlation Functions (HiCF and HDCF, respectively) have been implemented, following the work of M. Tuckerman et al. [[2]](https://aip.scitation.org/doi/10.1063/1.3474625) and T. Zelovich et al. [[3]](https://pubs.acs.org/doi/10.1021/acs.jpcc.8b10298).  
MD studies have revealed that the involved species (protons in this case) can be transferred back and forth between donors and acceptors in the femtoscale domain, and do not contribute to the overall transfer. This is known as "rattling" in the specialised literature [[2]](https://aip.scitation.org/doi/10.1063/1.3474625). To exclude this effect in the analysis, the standard practice is to remove snapshots from the trajectory. This can be a tedious and cumbersome task. To tackle this problem automatically, ALC_TRAJECTORY offers the ***rattling_wait*** directive, which handles the analysis as follows: once the transfer of the atomic species (in this case proton) from donor to acceptor has occurred, the algorithm waits/holds for 0.2 ps to accept the transfer. If the proton returns back to its donor, the transfer is discarded. If the ***rattling_wait*** directive is not declared, rattling effects will be included in the analysis.  
Activation of the ***&lifetime*** block generates the TCF and the TCF_AVG files. File TCF contains the computed quantities for all 10 ps segments, according to the specification of the ***&data_analysis*** block. File TCF_AVG, in contrast, reports the computed average (AVG), standard deviation (STD), as well as AVG+STD and AVG-STD. 

### <span style="color: maroon">Mean Square Displacement (MSD)</span>
As for the RDF, the computation of MSDs is a common feature in most available software for MD data analysis. The MSD is only computed for the monitored species, for which the definition of the ***&monitored_species*** is compulsory. For this example:

***&msd***  
&nbsp;&nbsp;&nbsp;select &nbsp;&nbsp;xy  
&nbsp;&nbsp;&nbsp;pbc_xyz&nbsp;&nbsp;   T&nbsp;&nbsp;  T&nbsp;&nbsp;  T  
***&end_msd***  

The ***select*** directive instructs the code to consider selected components of the atomic positions to compute the MSD. For those users that are not familiar with the formula to compute the MSD, we refer to the Supporting Information of Ref. [[4]](https://pubs.acs.org/doi/10.1021/acs.jpclett.1c04071?ref=PDF). In this example, the option "xy" computes the MSD only in the plane parallel to the Nafion membrane. Available options for ***select*** are: x, y, z, xy, xz, yz and xyz.  
As an optional directive, ***pbc_xyz*** allows including (or not) the effect of periodic boundary conditions (PBCs) for each Cartesian coordinate. PBCs are set by default.  
Activation of the ***&msd*** block generates the MSD and the MSD_AVG files. File MSD contains the computed quantities for all the 10 ps segments, according to the specification of the ***&data_analysis*** block. File MSD_AVG, in contrast, reports the computed average (AVG) and standard deviation(STD) as well as AVG+STD and AVG-STD. 

### <span style="color: maroon">Orientational Correlation Function (OCF)</span>
Together to the already presented capabilities, ALC_TRAJECTORY offers the computation of OCFs of monitored species, both for reactive and non-reactive systems. The monitored species must be a molecule (water in this case), otherwise the code will abort the execution. At time zero, the vector unit **u**(0) is computed using a chosen geometry criterion for each the monitored species (excluding the chemical species). Implemented options for geometry criterion are discussed briefly. Unit vectors define the orientations of the molecules. Likewise, at trajectory time **t** the unit vector **u**(t) is also computed for each monitored species. The OCT at time **t** is computed as `<`P<sub>l</sub> [**u**(t)**.u**(0)]`>`, where l is the order of the Legendre polynomial, the dot is the inner product and the average < > is done over the number of monitored species [[5]](https://pubs.acs.org/doi/10.1021/acs.jpcb.5b02936). A remarkable feature of the implementation is the possibility to compute OCF for reactive systems. The implemented algorithm identifies all the monitored species at time zero. If one of the species become a reactive site (from having accepted an atomic species), such site (formerly Ow and now Ow*) is discarded. Likewise, if the same site donates the species and becomes a monitored species again, it will be discarded for the whole time interval, as the correlation between **u**(t) and **u**(0) is already lost when the monitored species became an reactive site. In this example, we define the settings for OCF as follows:

***&ocf***  
&nbsp;&nbsp;&nbsp;legendre_order&nbsp;&nbsp;  2  
&nbsp;&nbsp;&nbsp;u_definition&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;    Plane  
***&end_ocf***  

The order of the polynomial "l" is set with the ***legendre_order*** directive. Although the OCF analysis is often performed by setting l=2, the user is free to select orders between 1 and 4.
To define the vector unit **u**, different geometrical criteria are implemented. The **u_definition** directive can be set to:

* bond_12: the unit vector between atoms 1 and 2 of the monitored species. This is the only option for diatomic molecules.
* bond_13: the unit vector between atoms 1 and 3 of the monitored species. 
* plane: the unit vector from the cross product between vector_12 and vector_13.
* bond_123: the unit vector from the **sum** of vector_12 and vector_13. 
* bond_12-13: the unit vector_12 and vector_13 are evaluated together within the average. 

It is users responsibility to test all these possible settings for the interpretation of the computed OCF. IMPORTANT: For monitored species with more than 2 atoms we do not recommend the use of options ***bond_12***, ***bond_13*** and ***plane***. 

Activation of the ***&ocf*** block generates the OCF and the OCF_AVG files. File OCF contains the computed quantities for all segments, according to the specification of the ***&data_analysis*** block. File OCF_AVG, in contrast, reports the computed average  (AVG) and standard deviation(STD) as well as AVG+STD and AVG-STD.

### <span style="color: maroon">Constraining the analysis to a selected region </span>
ALC_TRAJECTORY also offers the possibility to compute RDFs, MSDs, OCFs for a particular region of the simulation cell. If this block is omitted, the whole simulation cell is considered. To set a region, the user must set the ***&region*** block:

***&region***  
&nbsp;&nbsp;&nbsp;Delta_x &nbsp;&nbsp;&nbsp;&nbsp; -1.6 &nbsp;&nbsp; 11.6 &nbsp;&nbsp; inside  
&nbsp;&nbsp;&nbsp;Delta_y &nbsp;&nbsp;&nbsp;&nbsp;  -1.0 &nbsp;&nbsp; 14.0&nbsp;&nbsp;  inside  
&nbsp;&nbsp;&nbsp;Delta_z &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;   8.0 &nbsp;&nbsp; 12.0&nbsp;&nbsp;  outside  
***&end_region***  

The choice of the settings for ***Delta_x***, ***Delta_y*** and ***Delta_x*** is subject to the dimensions of the simulation cell. The first and second values are the minimum and maximum values for the domain along the corresponding coordinate in units of Angstrom. The third argument indicates if the region of interest for analysis are within (inside) or outside the given range. By comparison with the definition of the ***simulation_cell*** block above, we realise that ***Delta_x*** and ***Delta_y*** are redundant for this case, and they can be omitted. In fact, when the **Delta** setting is omitted, ALC_TRAJECTORY will use the whole spatial domain of the simulation cell for that coordinate. In this example, the analysis is focused on the interface regions (up and down) near the SO<sub>3</sub></sup>-</sup> groups. The central region between 8 and 12 Angstrom is not considered. If any of the defined settings is not consistent with the simulation cell (for example, replacing inside by outside in the specification of the ***Delta_x*** directive), the code will abort the execution. The user can double check the definition of the region in the generated OUTPUT file. Finally, multiple definitions for ***Delta_x***, ***Delta_y*** and ***Delta_z*** are allowed.

### <span style="color: maroon"> Tracking non-reactive species (unchanged chemistry) </span>
In addition to the tracking of the changing chemical species, the user can also choose to track atomic sites that do not change their chemistry. This can be important to compare how the location of reactive and non-reactive sites are distributed along the trajectory. The block must be defined as follows:

***&track_unchanged_chemistry***  
&nbsp;&nbsp;&nbsp; number &nbsp;&nbsp;&nbsp;4  
&nbsp;&nbsp;&nbsp; tag  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Sch  
&nbsp;&nbsp;&nbsp; list_indexes &nbsp;&nbsp;&nbsp; 133&nbsp;&nbsp;&nbsp;  134&nbsp;&nbsp;&nbsp; 135 &nbsp;&nbsp;&nbsp; 136  
***&end_track_unchanged_chemistry***  

The first directive must be ***number***, which indicates how many sites the user wants to track. A maximum of 10 sites are allowed. The directive ***tag*** specifies the atomic species to be tracked. Finally, the user must define the 4 atomic indexes with the ***indexes*** directive. In case that the declared indexes do not correspond to the defined ***tag***, the code will abort the execution. The positions of the tracked indexes are printed to the UNCHANGED_CHEMISTRY file.  
In the hypothetical scenario the users aims to track Ow atoms for this example, the code would track the Ow species as long as they do not change their chemistry. If the chemistry changes, results would be printed up to the frame where the  change has occurred.
