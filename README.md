# Optimization-algorithm-in-MATLAB

The codes that I wrote in my 4th academic year. (Intro to Optimization Class).

**Example from Steepest Descent Folder:**
Consider we need to find the local minimum x* of the objective function (In this case, Rosenbrock's function);
<img width="501" height="35" alt="image" src="https://github.com/user-attachments/assets/56bfdbf5-b61e-45cf-a713-eb77c0baaee8" />
By using modified Newton's Method to find the solution. The line search's search direction is <img width="221" height="31" alt="image" src="https://github.com/user-attachments/assets/92c7c7c5-6fe5-430e-8590-eb71143942fa" />. In case the Hessian matrix of the objective function is a singular matrix (the term <img width="72" height="39" alt="image" src="https://github.com/user-attachments/assets/3e583a1a-cdcf-4ce9-ab3b-878f6de2fca4" /> cannot be found.), We change the search direction method to steepest descent <img width="128" height="29" alt="image" src="https://github.com/user-attachments/assets/32f6c1b3-ca6d-4fc2-8157-937e549c583d" />.

Quick reminder: ? What is a line search and search direction?
<img width="1364" height="544" alt="image" src="https://github.com/user-attachments/assets/915c14e4-5729-407a-915e-2a292d48f94a" />

Here is the result. (After 51 iterations, we find the correct solution x = (1,1) with 10^-4 precision).
<img width="651" height="277" alt="image" src="https://github.com/user-attachments/assets/18cb2e33-c160-4d97-9786-9c192fbb5298" />

