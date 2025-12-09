This repository contains the codes and data to replicate the numerical results presented in the following CMAME paper.

Zhao J, Yoon H, Youn B D. An efficient concurrent topology optimization approach for frequency response problems[J]. Computer Methods in Applied Mechanics and Engineering, 2019, 347: 700-734.

[https://link.springer.com/article/10.1007/s00158-018-2044-x](https://www.sciencedirect.com/science/article/pii/S0045782519300088)

The purpose of this work is to develop an efficient concurrent topology optimization approach for minimizing frequency response of two-scale hierarchical structures over a given frequency interval. Compared with static problems, frequency response problems usually involve many load steps, which may lead to intensive computational burdens in both frequency response analysis and sensitivity analysis. This study thus proposes an enhanced decoupled sensitivity analysis method for frequency response problems, which is efficient even when plenty of frequency steps are involved and/or damping is considered. Furthermore, a combined method of modal superposition and model order reduction is incorporated for efficient frequency response analysis of two-scale hierarchical structures. A modified threshold Heaviside projection method is used to obtain black-and-white designs and the method of moving asymptotes (MMA) is employed to update the design variables. Several numerical examples are presented to demonstrate the effectiveness of the proposed approach.

Fig. 1. Two-scale concurrent topology optimization framework with homogeneous porous material

![1-s2 0-S0045782519300088-gr1_lrg](https://github.com/user-attachments/assets/4648c234-21df-420a-8557-b28c6d68030f)

Fig. 2. Problem definition of concurrent topology optimization for minimizing displacement response.

![1-s2 0-S0045782519300088-gr3_lrg](https://github.com/user-attachments/assets/f21c6fc7-7617-44ac-80c6-fe2bceb35037)

Fig. 3. Comparison of sensitivities obtained by the original and enhanced decoupled methods.

![1-s2 0-S0045782519300088-gr4_lrg](https://github.com/user-attachments/assets/21b0eaff-da0e-4794-be55-06f4c398f3df)

Fig. 4. CPU times of frequency response analysis and the two decoupled sensitivity analysis methods.

![1-s2 0-S0045782519300088-gr5_lrg](https://github.com/user-attachments/assets/85516294-80fc-4d19-98fb-7b95f53bbb92)

Fig. 5. Concurrent optimized designs of the cantilever beam for minimizing displacement response, [0,10π]Hz.

![1-s2 0-S0045782519300088-gr7_lrg](https://github.com/user-attachments/assets/5786504e-f94a-4b8d-96a9-179f5f558e96)

Fig. 6. Concurrent optimized designs of the cantilever beam for minimizing displacement response, [0,100π]Hz.
![1-s2 0-S0045782519300088-gr8_lrg](https://github.com/user-attachments/assets/ad9f3c19-47e4-470a-9a49-81114a8cf517)

Fig. 7. Concurrent optimized designs of the cantilever beam for minimizing displacement response, [0,150π]Hz.
![1-s2 0-S0045782519300088-gr9_lrg](https://github.com/user-attachments/assets/91dd22a7-664d-4951-b359-dbcf7a395174)

Fig. 8. Displacement response curves for the cantilever beam designs.
![1-s2 0-S0045782519300088-gr13_lrg](https://github.com/user-attachments/assets/aa6fcf32-c10d-471c-8096-7c3b5b674ad5)

Fig. 9. Problem definition of concurrent topology optimization for minimizing acceleration response.
![1-s2 0-S0045782519300088-gr15_lrg](https://github.com/user-attachments/assets/5c5f9a47-8887-4b4f-b933-ee578414f2e5)

Fig. 10. Concurrent optimized designs of the cantilever beam for minimizing acceleration response, [0,10π]Hz.
![1-s2 0-S0045782519300088-gr16_lrg](https://github.com/user-attachments/assets/2ca98896-d2a0-456e-8a58-c82aff87b071)

Fig. 11. Concurrent optimized designs of the cantilever beam for minimizing acceleration response, [0,100π]Hz.
![1-s2 0-S0045782519300088-gr17_lrg](https://github.com/user-attachments/assets/740973ff-2765-4ce6-a34a-16533bb75084)

Fig. 12. Concurrent optimized designs of the cantilever beam for minimizing acceleration response, [0,150π]Hz.
![1-s2 0-S0045782519300088-gr18_lrg](https://github.com/user-attachments/assets/726505d8-640a-4267-a339-2174c61a64e2)

Fig. 13. Acceleration response curves for the cantilever beam designs.
![1-s2 0-S0045782519300088-gr22_lrg](https://github.com/user-attachments/assets/fabe375e-ae46-40a9-b741-4021b64721cd) 

Fig. 14. Concurrent optimized designs of the cantilever beam for minimizing acceleration response without constraint on the static stiffness, [0,10π]Hz.
![1-s2 0-S0045782519300088-gr23_lrg](https://github.com/user-attachments/assets/586fca8f-9180-457b-92ce-30f2b561a255)

Fig. 15. Iteration history of concurrent design of the cantilever beam for minimizing acceleration response without constraint on the static stiffness, [0,100π]Hz. 
![1-s2 0-S0045782519300088-gr24_lrg](https://github.com/user-attachments/assets/1523584c-92ca-4f5d-9dd8-9a256f71a05d)

