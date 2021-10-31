Data File Format of rcvr.dat:

rcvr.dat is an 8x7 matrix containing raw ranging information. Each of the 8 rows contains independent measurements for each of 8 satellites in view at the current epoch (an epoch is simply a term refers to a single discrete time; since our receivers provide data at approximately 1 sec. intervals, each epoch occurs approximately 1 sec. after the prior epoch. The columns of this matrix in clued the following data:
<img width="718" alt="Screenshot 2021-10-31 at 12 53 22 PM" src="https://user-images.githubusercontent.com/71690213/139568405-2df2e8bc-b392-466f-a203-a5cfef60ffce.png">

Data File Format of eph.dat:

eph.dat is an 8 x 24 matrix containing the ephemeris data from a GPS receiver. This data is used to estimate the orbital position of each satellite at any given time. Each row contains ephemeris data for a single satellite. The columns of this matrix include the following data:

<img width="718" alt="Screenshot 2021-10-31 at 12 54 45 PM" src="https://user-images.githubusercontent.com/71690213/139568436-749bfd48-e6e9-4159-88d5-a9e850e71764.png">

