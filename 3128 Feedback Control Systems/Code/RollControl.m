function aileron_perturb = RollControl(ka, kp,phi_c,phi,p)
    aileron_perturb = ka*(kp*(phi_c - phi) - p); 
end