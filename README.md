# Free Vibration of Laminated Composite Plates

Hi. This repository contains MATLAB codes for analyzing the free vibration of laminated composite plates. 

## Sources & Inspiration

The foundation of this work is based on Ferreira's MATLAB Codes for Finite Element Analysis [1]. The code in the book is simple, but as a beginner, I found the descriptions a bit confusing.

To better understand the engineering formulation (especially for the stiffness K and mass M matrices), I relied on the works of Pervez and Zabaras [2][3]. Even though their work focuses on damping and calculating C, their formulation for K and M is excellent. Combining their theory with Ferreira's code really helped me grasp the model fully.

## The Code (and Modifications)

As mentioned, the core of this project is based on Ferreira's book, so don't be surprised if you see a lot of copy-pasting!😅

However, I did add some missing parts and modified the code to make it more suitable for my own understanding and future projects.

Key Modification:
I changed the Degree of Freedom (DoF) ordering.

* Original: Grouped by variable (e.g., w1, w2..., then fix1, fix2..., etc.)
* My version: Node-wise ordering (e.g., u1 v1 w1 fix1 fiy1, then u2 v2 w2...)

## Known Issues (Help!)
I'm facing a specific issue that I haven't been able to solve yet:

* Singularity in K Matrix: When I use 8 or 9-node elements, I get a singularity in my Stiffness (K) matrix.
* What I tried: I tried adding a proportion of the Mass matrix to M and K, but it didn't fix the issue.

If you know why this happens or how to fix it, I would appreciate your help! Also, if you find other bugs that I don't even know about, fix it...(but don't tell anyone😅)

## References 
[1] Ferreira, António JM. MATLAB codes for finite element analysis: solids and structures. Dordrecht: Springer Netherlands, 2009.

[2] Pervez, Tasneem, and Nicholas Zabaras. "Transient dynamic and damping analysis of laminated anisotropic plates using a refined plate theory." International journal for numerical methods in engineering 33.5 (1992): 1059-1080.

[3] Zabaras, Nicholas, and Tasneem Pervez. "Viscous damping approximation of laminated anisotropic composite plates using the finite element method." Computer methods in applied mechanics and engineering 81.3 (1990): 291-316.
