# QE band structure plotting tools
Does not require running bands.x within QE
Directly plot band structure after 'nscf'/'bands' calculation 

'kpath' file needs to be written containing Kpoints  in the following format
All high-symmetry kpoints in the band structure should be in 'kpath' file  
Label   Kx  Ky Kz
e.g.
X     0.5   0.0  0.0

