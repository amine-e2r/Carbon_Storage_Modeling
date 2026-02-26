## Project Description  

Numerical modeling of carbon dioxide storage using a system of nonlinear differential equations to simulate carbon exchanges between the atmosphere, trees, and soil compartments.  

- Developed a coupled dynamical system describing carbon fluxes between atmosphere ($C_A$), vegetation ($C_T$), and soil ($C_S$).  
- Modeled sequestration through a logistic growth term controlled by parameters $\alpha$ (sequestration rate) and $K$ (maximum storage capacity).  
- Incorporated biological transfer mechanisms:  
  - $\beta$ (tree respiration to atmosphere)  
  - $\gamma$ (litter transfer to soil)  
  - $\delta$ (soil respiration and tree-to-soil transfer)  
- Implemented and compared three numerical schemes:  
  - Explicit Euler  
  - Implicit Euler (fixed-point resolution)  
  - Trapezoidal method with Newton iteration  
- Evaluated numerical stability, convergence behavior, and sensitivity to time step size.  
- Extended the model by adding an ocean compartment and an anthropogenic emission term to improve realism.  

---

## Key Results  

- **Trapezoidal method with Newton** demonstrated the highest numerical stability and remained reliable for large time steps.  
- **Implicit Euler** improved stability compared to the explicit scheme but required iterative resolution due to nonlinearity.  
- **Explicit Euler** became unstable when the time step was too large, particularly during rapid dynamic transitions.  
- Sensitivity analysis showed that **$\alpha$ (sequestration rate)** is the most influential parameter for reducing atmospheric CO₂.  
- Higher values of **$\beta$ and $\delta$** increased atmospheric carbon levels and could prevent the system from reaching equilibrium.  
- The long-term behavior of the system strongly depends on parameter interactions, leading either to stable equilibrium, enhanced sequestration, or atmospheric carbon growth.  

This project highlights how nonlinear dynamical systems combined with appropriate numerical methods can effectively simulate forest carbon storage dynamics and emphasizes the importance of numerical stability in environmental modeling.
