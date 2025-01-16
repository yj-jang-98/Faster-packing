# CDSL: Cryptography for Dynamic Systems Library

Still developing readme section

CDSL provides codes for implementing secure dynamic systems based on modern cryptography.
The library features linear dynamic controllers operating over homomorphically encrypted data implemented using [Lattigo](https://github.com/tuneinsight/lattigo) version 6.1.0.
The encrypted controllers are designed based on the state-of-the-art methods developed by CDSL, [SNU](https://post.cdsl.kr/) and [SEOULTECH](https://junsookim4.wordpress.com/).


---

### Overview


Given a plant 

$$
\begin{aligned}
x_p(t+1) &= Ax_p(t) + Bu(t), \quad x_p(0) = x_p^{\mathsf{ini}} \\
y(t) &= Cx_p(t)
\end{aligned}
$$

and a pre-designed stabilizng controller (which is controllable and observable)

$$
\begin{aligned}
x(t+1) &= Fx(t) + Gy(t), \quad x(0) = x^{\mathsf{ini}} \\
u(t) &= Hx(t)
\end{aligned}
$$

this code provides two methods to operate the linear dynamic controller over encrypted data, using a Ring-LWE based cryptosystem. 


- `ctrRGSW` [1]: Supports unlimited number of recursive homomorhpic multiplications without the use of bootstrapping. More specifically, the encrypted controller state is recursively multiplied to the encrypted state matrix without decryption. The effect of error growth is suppressed by the stability of the closed-loop system. 
    - `ctrRGSW/noPacking`: Naive implementation that does not use packing. 
    - `ctrRGSW/packing`: A novel "coefficient packing" technique is applied, resulting in enhanced computation speed and memory efficiency   
    - `ctrRGSW/conversion.m`: Converts the state matrix of the controller into integers based on the apporach of [3]:
       - Given $F$ and $H$, it finds an appropriate $R$ such that $F-RH$ is an integer matrix. Then, the state dynamics of the controller can be rewritten as follows regarding $u(t)$ as a fed-back input.
       - One may use other methods, such as [4-6], to convert the state matrix into integers while preserving the control performance without using re-encryption. 


$$
x(t+1) = (F-RH)x(t) + Gy(t) + Ru(t)
$$




- `ctrRLWE` [2]: 

- In this example, we use the four-tank system [7] as the plant and pre-design an observer-based controller. Please refer to [8] for more details regarding the setting.

---

### How to use
Download or clone this repository using
```
git clone https://github.com/CDSL-EncryptedControl/CDSL.git
```


Then, change the directory to a folder you wish to use and run the `main.go` file. For example,

```
cd ctrRGSW/noPacking
go run main.go  
```
or
```
cd ctrRGSW/packing
go run main.go  
```
or
```
cd ctrRLWE
go run main.go  
```
on the terminal.


---

### Set parameters 

* `rlwe.NewParametersFromLiteral`: Ring-LWE parameters (LogN = 11 and LogQ = 54 gives $N=2^{11}$ and some prime $q$ such that $q \approx 2^{54}$)

* `s`, `L`, and `r`: Scale factors 

* `iter`: Number of iterations for simulation 

* `A`, `B`, and `C`: State space matrices of the discrete time plant written by

> $x(t+1) = Ax(t) + Bu(t), \quad y(t) = Cx(t)$

* `F`, `G`, `R` and `H`: State space matrices of the discrete time controller. 
Given a controller of the form 
> $x(t+1) = Kx(t) + Gy(t), \quad u(t) = Hx(t)$

one can regard $u(t)$ as a fed-back input and design $R$, so that the state matrix $F:=K-RH$ of
> $x(t+1) = (K-RH)x(t) + Gy(t)+Ru(t), \quad u(t) = Hx(t)$

consists of integers. More details can be found in Section 5 of [1] or Lemma 1 of [2].

* `xPlantInit`, `xContInit`: Initial conditions of the plant and the controller

* `tau`: Least power of two greater than the dimensions of the state, output, and input of the controller (Only used in `Ring-GSW_Packed.go`)

---

### References
[1] [Y. Jang, J. Lee, S. Min, H. Kwak, J. Kim, and Y. Song, "Ring-LWE based encrypted controller with unlimited number of recursive multiplications and effect of error growth," 2024, arXiv:2406.14372.](https://arxiv.org/abs/2406.14372)

[2] [J. Lee, D. Lee, J. Kim, and H. Shim, "Encrypted dynamic control exploiting limited number of multiplications and a method using RLWE-based cryptosystem," _IEEE Trans. Syst. Man. Cybern.: Syst._, vol. 55, no. 1, pp. 158-169, 2025.](https://ieeexplore.ieee.org/abstract/document/10730788)

[3] [J. Kim, H. Shim, and K. Han, "Dynamic controller that operates over homomorphically encrypted data for infinite time horizon," _IEEE Trans. Autom. Control_, vol. 68, no. 2, pp. 660-672, 2023.](https://ieeexplore.ieee.org/abstract/document/9678042)

[4] [J. Kim, H. Shim, H. Sandberg, and K. H. Johansson, “Method for running dynamic systems over encrypted data for infinite time horizon without bootstrapping and re-encryption,” in Proc. 60th IEEE Conf. Decision Control, 2021, pp. 5614–5619.](https://ieeexplore.ieee.org/abstract/document/9682828?casa_token=LHR79rToQ7oAAAAA:Wz1AzFWR7VW6DYKUhLFYcoXtpMx4AIT9E_krpOpFy7QUO5lSkvPf_0ZZgPsdp65ZzaGx-ejlPA)

[5] [M. S. Tavazoei, “Non-minimality of the realizations and possessing state matrices with integer elements in linear discrete-time controllers,” IEEE Trans. Autom. Control, vol. 68, no. 6, pp. 3698–3703, 2023.](https://ieeexplore.ieee.org/abstract/document/9835020?casa_token=_rdGjQLc7ZEAAAAA:QLxzC1QlnNVYriMTL1gbSjtv5U2oTwfVO5OqVFfGS0Qpz8hx7exSuJKJ9H8XBh_qDucoZt8oBg)

[6] [J. Lee, D. Lee, S. Lee, J. Kim, and H. Shim, “Conversion of controllers to have integer state matrix for encrypted control: Non-minimal order approach,” in Proc. 62nd IEEE Conf. Decision Control, 2023, pp. 5091–5096.](https://ieeexplore.ieee.org/abstract/document/10383200?casa_token=lbob37tAZ-MAAAAA:vAVUmuIngRzHYefqaYHQM5TfukcAI7Lh1YmYngqcLYMj74Mtzq0xGybkntfWSd-DKwogxrvnxg)

[7] [K. Johansson, “The quadruple-tank process: A multivariable laboratory process with an adjustable zero,” IEEE Trans. Control Sys. Technol., vol. 8, no. 3, pp. 456–465, 2000.](https://ieeexplore.ieee.org/abstract/document/845876?casa_token=1CWEIgmKIscAAAAA:Hh3D4_xn5B8MWVoMpQHof8glwtWpGXMuddehBoKXbZAOh2WwsDlemeiWeZ6nAwQGThjhYYw1wQ)

[8] [Y. Jang, J. Lee, and J. Kim, "Documentation on encrypted dynamic control simulation code using Ring-LWE based cryptosystems," ](link)
