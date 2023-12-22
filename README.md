**Code used for analysis for [Spatial variation in drainage area —
    Runoff relationships and implications for bankfull geometry scaling](https://www.sciencedirect.com/science/article/pii/S0169555X2300418X?casa_token=Lxa1OVNejeYAAAAA:2RlnyR2jTHuEgApmhp_0eBk4lYklhdAKPf4a8q_Prr2vfiDxH-T2Eq0TF7jUAu-TlelEYE7SZw)**

*Citation:*

Erikson, C. M., Renshaw, C. E., & Magilligan, F. J. (2024). Spatial variation in drainage area —
    Runoff relationships and implications for bankfull geometry scaling. *Geomorphology*, 446,
    108998. https://doi.org/10.1016/j.geomorph.2023.108998

# Motivation
River discharge is a fundamental component of a wide range of analyses, including flood frequency assesments and landscape evolution. However, discharge observations are both spatially and temporally confined by the existence of stream gauges.
To obviate this issue, discharge is calcualted as a power law of watershed drainage area. 

<img width="100" alt="image" src="https://github.com/cmerikson/Runoff_Scaling/assets/109803481/6bb31b1b-72ef-4945-bc47-db5ea3f2baf7">

A common assumption in this calculation is that drainage area and discharge have a linear relationship (c=1). 
We test that assumption across the United states and Canada. We find that, while most of the continent does indeed demonstate a linear relationship between discharge and drainage area, places like the American Mid-West and Southwest have non-linear relationships.
We also find that locations with non-linear scaling tend to have low runoff efficiency. Stream channels in these watersheds may be primarily shaped by larger flows than channels in other watersheds.

# Discharge - Drainage Area Scaling
![A kriged map of the exponent describing the scaling of drainage area with discahrge.](https://github.com/cmerikson/Runoff_Scaling/files/13747156/Q2_map.pdf)
