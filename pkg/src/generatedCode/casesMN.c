#define psiexp(i, j) psiexp_array[3*(i - -3*q) + j - 1]
#define Psi(j) psiParam[j - 1]
const scalar h = *hParam;
scalar *psiexp_array = malloc(sizeof(scalar)*3*(1 + 3*q + Max(1 + n - q,4 + 2*q)));

for (int i = -3*q; i <= Max(1 + n - q,4 + 2*q); i++)
{
	for (int j = 1; j <= 3; j++)
	{
		psiexp(i, j) = Exponential(h*i*Psi(j));
	}
}

for (int j = 0; j <= 0; j++)
{
	for (int k = 0; k <= 0; k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = 0; nu <= -2 + q; nu++)
			{
				scalar sumInuP1;
				{
					scalar sumIkappa = 0;
					for (int kappa = Max(0,2 + nu); kappa <= q; kappa++)
					{
						sumIkappa += (vEG2G2c*(2*n - 2*q - Min(kappa - nu,n - q))*(-1 + Min(kappa - nu,n - q))*b[kappa])/2. + vEG2G2cAW*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(kappa - nu - Min(kappa - nu,n - q),1)*((1 + n - q)*psiexp(1,1) + Min(kappa - nu,n - q)*(-1 + psiexp(1,1))*psiexp(1,1) + (-n + q)*psiexp(2,1) + (-n + q)*psiexp(Min(kappa - nu,n - q),1) + (-1 + n - q)*psiexp(1 + Min(kappa - nu,n - q),1));
					}
					sumInuP1 = sumIkappa;
				}
				sumInu += sumInuP1*b[nu];
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += MMjkP1/(n - q);
	}
}

for (int j = 0; j <= 0; j++)
{
	for (int k = 0; k <= 0; k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = 0; nu <= q; nu++)
			{
				scalar sumInuP1;
				{
					scalar sumIkappa = 0;
					for (int kappa = 0; kappa <= Min(-2 + n + nu - q,q); kappa++)
					{
						sumIkappa += (vEG2G2c*b[kappa]*(-n + q - 2*n*q + (1 - 2*n + 2*q)*Max(0,kappa - nu) + Power2(n) + Power2(q) + Power2(Max(0,kappa - nu))))/2. + vEG2G2cAW*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-kappa + nu,1)*(psiexp(1 + n - q,1) + (-1 + n - q)*psiexp(1 + Max(0,kappa - nu),1) + Max(0,kappa - nu)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,kappa - nu),1) + (-n + q)*psiexp(2 + Max(0,kappa - nu),1));
					}
					sumInuP1 = sumIkappa;
				}
				sumInu += sumInuP1*b[nu];
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += MMjkP1/(n - q);
	}
}

for (int j = 0; j <= 0; j++)
{
	for (int k = 0; k <= 0; k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = 0; nu <= -1 + q; nu++)
			{
				scalar sumInuP1;
				{
					scalar sumIkappa = 0;
					for (int kappa = Max(0,1 + nu); kappa <= Min(-1 + n + nu - q,q); kappa++)
					{
						sumIkappa += (-kappa + n + nu - q)*b[kappa];
					}
					sumInuP1 = sumIkappa;
				}
				sumInu += sumInuP1*b[nu];
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += (MMjkP1*vEG4c)/(n - q);
	}
}

for (int j = 0; j <= 0; j++)
{
	for (int k = 0; k <= 0; k++)
	{
		scalar MMjkP1;
		{
			scalar sumIkappa = 0;
			for (int kappa = 0; kappa <= q; kappa++)
			{
				sumIkappa += b[kappa];
			}
			MMjkP1 = sumIkappa;
		}
		MM(j, k) += (azero*MMjkP1*vEG2c*(-n + q - 2*n*q + Power2(n) + Power2(q)))/(n - q);
	}
}

for (int j = 0; j <= 0; j++)
{
	for (int k = 0; k <= 0; k++)
	{
		MM(j, k) += ((-1 + n - q)*Power2(azero))/2.;
	}
}

for (int j = 1; j <= -3 + q; j++)
{
	for (int k = 0; k <= 0; k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = 1 + j; nu <= -2 + q; nu++)
			{
				scalar sumInuP1;
				{
					scalar sumIkappa = 0;
					for (int kappa = Max(0,2 + nu); kappa <= q; kappa++)
					{
						sumIkappa += (vEG2G2G2c*(2*n - 2*q - Min(kappa - nu,n - q))*(-1 + Min(kappa - nu,n - q))*b[kappa])/2. + (vEG2G2G2cAW*(2*n - 2*q - Min(kappa - nu,n - q))*(-1 + Min(kappa - nu,n - q))*b[kappa]*psiexp(-j + nu,1))/2. + vEG2G2G2cAWAD*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-j + kappa - Min(kappa - nu,n - q),1)*((1 + n - q)*psiexp(1,1) + Min(kappa - nu,n - q)*(-1 + psiexp(1,1))*psiexp(1,1) + (-n + q)*psiexp(2,1) + (-n + q)*psiexp(Min(kappa - nu,n - q),1) + (-1 + n - q)*psiexp(1 + Min(kappa - nu,n - q),1)) + vEG2G2G2cAD*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(kappa - nu - Min(kappa - nu,n - q),1)*((1 + n - q)*psiexp(1,1) + Min(kappa - nu,n - q)*(-1 + psiexp(1,1))*psiexp(1,1) + (-n + q)*psiexp(2,1) + (-n + q)*psiexp(Min(kappa - nu,n - q),1) + (-1 + n - q)*psiexp(1 + Min(kappa - nu,n - q),1)) + vEG2G2G2cAWBD*b[kappa]*Power2(1/(-1 + psiexp(1,2)))*psiexp(-j + nu,1)*psiexp(kappa - nu - Min(kappa - nu,n - q),2)*((1 + n - q)*psiexp(1,2) + Min(kappa - nu,n - q)*(-1 + psiexp(1,2))*psiexp(1,2) + (-n + q)*psiexp(2,2) + (-n + q)*psiexp(Min(kappa - nu,n - q),2) + (-1 + n - q)*psiexp(1 + Min(kappa - nu,n - q),2));
					}
					sumInuP1 = sumIkappa;
				}
				sumInu += sumInuP1*b[nu];
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += MMjkP1/(n - q);
	}
}

for (int j = 1; j <= -2 + q; j++)
{
	for (int k = 0; k <= 0; k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = Max(2 + j,2 - n + q); nu <= q; nu++)
			{
				scalar sumInuP1;
				{
					scalar sumIkappa = 0;
					for (int kappa = Max(0,2 + j); kappa <= Min(-2 + n + nu - q,q); kappa++)
					{
						sumIkappa += -(vEG2G2G2c*(2*n - 2*q - Max(0,kappa - nu) - Min(-j + kappa,n - q))*(1 + Max(0,kappa - nu) - Min(-j + kappa,n - q))*b[kappa])/2. - (vEG2G2G2cAWAD*(2*n - 2*q - Max(0,kappa - nu) - Min(-j + kappa,n - q))*(1 + Max(0,kappa - nu) - Min(-j + kappa,n - q))*b[kappa]*psiexp(-j + nu,1))/2. + vEG2G2G2cAW*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-j + kappa - Max(0,kappa - nu) - Min(-j + kappa,n - q),1)*(psiexp(1 + Max(0,kappa - nu),1) + n*psiexp(1 + Max(0,kappa - nu),1) - q*psiexp(1 + Max(0,kappa - nu),1) + Min(-j + kappa,n - q)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,kappa - nu),1) - n*psiexp(2 + Max(0,kappa - nu),1) + q*psiexp(2 + Max(0,kappa - nu),1) - n*psiexp(Min(-j + kappa,n - q),1) + q*psiexp(Min(-j + kappa,n - q),1) - Max(0,kappa - nu)*(-1 + psiexp(1,1))*psiexp(Min(-j + kappa,n - q),1) - psiexp(1 + Min(-j + kappa,n - q),1) + n*psiexp(1 + Min(-j + kappa,n - q),1) - q*psiexp(1 + Min(-j + kappa,n - q),1)) + vEG2G2G2cAD*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-kappa + nu,1)*(-psiexp(1 + Max(0,kappa - nu),1) + n*psiexp(1 + Max(0,kappa - nu),1) - q*psiexp(1 + Max(0,kappa - nu),1) + Max(0,kappa - nu)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,kappa - nu),1) - n*psiexp(2 + Max(0,kappa - nu),1) + q*psiexp(2 + Max(0,kappa - nu),1) - n*psiexp(Min(-j + kappa,n - q),1) + q*psiexp(Min(-j + kappa,n - q),1) - Min(-j + kappa,n - q)*(-1 + psiexp(1,1))*psiexp(Min(-j + kappa,n - q),1) + psiexp(1 + Min(-j + kappa,n - q),1) + n*psiexp(1 + Min(-j + kappa,n - q),1) - q*psiexp(1 + Min(-j + kappa,n - q),1)) + vEG2G2G2cAWBD*b[kappa]*Power2(1/(psiexp(1,1) - psiexp(1,2)))*psiexp(-kappa + nu,2)*psiexp(-j + kappa - Max(0,kappa - nu) - Min(-j + kappa,n - q),1)*(Max(0,kappa - nu)*(-psiexp(1,1) + psiexp(1,2))*psiexp(1 + Max(0,kappa - nu),2)*psiexp(Min(-j + kappa,n - q),1) - n*psiexp(2 + Max(0,kappa - nu),2)*psiexp(Min(-j + kappa,n - q),1) + q*psiexp(2 + Max(0,kappa - nu),2)*psiexp(Min(-j + kappa,n - q),1) + Min(-j + kappa,n - q)*(psiexp(1,1) - psiexp(1,2))*psiexp(1 + Max(0,kappa - nu),1)*psiexp(Min(-j + kappa,n - q),2) - n*psiexp(2 + Max(0,kappa - nu),1)*psiexp(Min(-j + kappa,n - q),2) + q*psiexp(2 + Max(0,kappa - nu),1)*psiexp(Min(-j + kappa,n - q),2) - psiexp(1 + Max(0,kappa - nu),2)*psiexp(1 + Min(-j + kappa,n - q),1) + n*psiexp(1 + Max(0,kappa - nu),2)*psiexp(1 + Min(-j + kappa,n - q),1) - q*psiexp(1 + Max(0,kappa - nu),2)*psiexp(1 + Min(-j + kappa,n - q),1) + psiexp(1 + Max(0,kappa - nu),1)*psiexp(1 + Min(-j + kappa,n - q),2) + n*psiexp(1 + Max(0,kappa - nu),1)*psiexp(1 + Min(-j + kappa,n - q),2) - q*psiexp(1 + Max(0,kappa - nu),1)*psiexp(1 + Min(-j + kappa,n - q),2));
					}
					sumInuP1 = sumIkappa;
				}
				sumInu += sumInuP1*b[nu];
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += MMjkP1/(n - q);
	}
}

for (int j = 1; j <= -1 + q; j++)
{
	for (int k = 0; k <= 0; k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = 1 + j; nu <= q; nu++)
			{
				scalar sumInuP1;
				{
					scalar sumIkappa = 0;
					for (int kappa = 0; kappa <= Min(-2 + j + n - q,q); kappa++)
					{
						sumIkappa += (vEG2G2G2c*b[kappa]*(-n + q - 2*n*q + (1 - 2*n + 2*q)*Max(0,-j + kappa) + Power2(n) + Power2(q) + Power2(Max(0,-j + kappa))))/2. + (vEG2G2G2cAD*b[kappa]*(-n + q - 2*n*q + (1 - 2*n + 2*q)*Max(0,-j + kappa) + Power2(n) + Power2(q) + Power2(Max(0,-j + kappa)))*psiexp(-j + nu,1))/2. + vEG2G2G2cAW*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(j - kappa,1)*(psiexp(1 + n - q,1) + (-1 + n - q)*psiexp(1 + Max(0,-j + kappa),1) + Max(0,-j + kappa)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,-j + kappa),1) + (-n + q)*psiexp(2 + Max(0,-j + kappa),1)) + vEG2G2G2cAWBD*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(j - kappa,1)*psiexp(-j + nu,2)*(psiexp(1 + n - q,1) + (-1 + n - q)*psiexp(1 + Max(0,-j + kappa),1) + Max(0,-j + kappa)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,-j + kappa),1) + (-n + q)*psiexp(2 + Max(0,-j + kappa),1)) + vEG2G2G2cAWAD*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-kappa + nu,1)*(psiexp(1 + n - q,1) + (-1 + n - q)*psiexp(1 + Max(0,-j + kappa),1) + Max(0,-j + kappa)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,-j + kappa),1) + (-n + q)*psiexp(2 + Max(0,-j + kappa),1));
					}
					sumInuP1 = sumIkappa;
				}
				sumInu += sumInuP1*b[nu];
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += MMjkP1/(n - q);
	}
}

for (int j = 1; j <= -2 + q; j++)
{
	for (int k = 0; k <= 0; k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = Max(1 + j,1 - n + q); nu <= -1 + q; nu++)
			{
				scalar sumInuP1;
				{
					scalar sumIkappa = 0;
					for (int kappa = Max(0,1 + nu); kappa <= Min(-1 + n + nu - q,q); kappa++)
					{
						sumIkappa += (-kappa + n + nu - q)*b[kappa];
					}
					sumInuP1 = sumIkappa;
				}
				sumInu += sumInuP1*b[nu]*(vEG2G4c + vEG2G4cAW*psiexp(-j + nu,1));
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += MMjkP1/(n - q);
	}
}

for (int j = 1; j <= -1 + q; j++)
{
	for (int k = 0; k <= 0; k++)
	{
		scalar MMjkP1;
		{
			scalar sumIkappa = 0;
			for (int kappa = Max(0,1 + j); kappa <= Min(-1 + j + n - q,q); kappa++)
			{
				sumIkappa += (j - kappa + n - q)*b[kappa];
			}
			MMjkP1 = sumIkappa;
		}
		scalar MMjkP2;
		{
			scalar sumInu = 0;
			for (int nu = 1 + j; nu <= q; nu++)
			{
				sumInu += vEG4G2c*b[nu] + vEG4G2cAD*b[nu]*psiexp(-j + nu,1) + vEG4G2cBD*b[nu]*psiexp(-j + nu,2);
			}
			MMjkP2 = sumInu;
		}
		MM(j, k) += (MMjkP1*MMjkP2)/(n - q);
	}
}

for (int j = 1; j <= -2 + q; j++)
{
	for (int k = 0; k <= 0; k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = 0; nu <= -1 + j; nu++)
			{
				scalar sumInuP1;
				{
					scalar sumIkappa = 0;
					for (int kappa = Max(0,2 + j); kappa <= q; kappa++)
					{
						sumIkappa += (vEG2G2G2c*(2*n - 2*q - Min(-j + kappa,n - q))*(-1 + Min(-j + kappa,n - q))*b[kappa])/2. + (vEG2G2G2cAW*(2*n - 2*q - Min(-j + kappa,n - q))*(-1 + Min(-j + kappa,n - q))*b[kappa]*psiexp(j - nu,1))/2. + vEG2G2G2cAD*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-j + kappa - Min(-j + kappa,n - q),1)*((1 + n - q)*psiexp(1,1) + Min(-j + kappa,n - q)*(-1 + psiexp(1,1))*psiexp(1,1) + (-n + q)*psiexp(2,1) + (-n + q)*psiexp(Min(-j + kappa,n - q),1) + (-1 + n - q)*psiexp(1 + Min(-j + kappa,n - q),1)) + vEG2G2G2cAWAD*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(kappa - nu - Min(-j + kappa,n - q),1)*((1 + n - q)*psiexp(1,1) + Min(-j + kappa,n - q)*(-1 + psiexp(1,1))*psiexp(1,1) + (-n + q)*psiexp(2,1) + (-n + q)*psiexp(Min(-j + kappa,n - q),1) + (-1 + n - q)*psiexp(1 + Min(-j + kappa,n - q),1)) + vEG2G2G2cAWBD*b[kappa]*Power2(1/(-1 + psiexp(1,2)))*psiexp(j - nu,1)*psiexp(-j + kappa - Min(-j + kappa,n - q),2)*((1 + n - q)*psiexp(1,2) + Min(-j + kappa,n - q)*(-1 + psiexp(1,2))*psiexp(1,2) + (-n + q)*psiexp(2,2) + (-n + q)*psiexp(Min(-j + kappa,n - q),2) + (-1 + n - q)*psiexp(1 + Min(-j + kappa,n - q),2));
					}
					sumInuP1 = sumIkappa;
				}
				sumInu += sumInuP1*b[nu];
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += MMjkP1/(n - q);
	}
}

for (int j = 2; j <= q; j++)
{
	for (int k = 0; k <= 0; k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = 0; nu <= Min(-2 + j,-2 + q); nu++)
			{
				scalar sumInuP1;
				{
					scalar sumIkappa = 0;
					for (int kappa = Max(0,2 + nu); kappa <= Min(-2 + j + n - q,q); kappa++)
					{
						sumIkappa += -(vEG2G2G2c*(2*n - 2*q - Max(0,-j + kappa) - Min(kappa - nu,n - q))*(1 + Max(0,-j + kappa) - Min(kappa - nu,n - q))*b[kappa])/2. - (vEG2G2G2cAWAD*(2*n - 2*q - Max(0,-j + kappa) - Min(kappa - nu,n - q))*(1 + Max(0,-j + kappa) - Min(kappa - nu,n - q))*b[kappa]*psiexp(j - nu,1))/2. + vEG2G2G2cAW*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(kappa - nu - Max(0,-j + kappa) - Min(kappa - nu,n - q),1)*(psiexp(1 + Max(0,-j + kappa),1) + n*psiexp(1 + Max(0,-j + kappa),1) - q*psiexp(1 + Max(0,-j + kappa),1) + Min(kappa - nu,n - q)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,-j + kappa),1) - n*psiexp(2 + Max(0,-j + kappa),1) + q*psiexp(2 + Max(0,-j + kappa),1) - n*psiexp(Min(kappa - nu,n - q),1) + q*psiexp(Min(kappa - nu,n - q),1) - Max(0,-j + kappa)*(-1 + psiexp(1,1))*psiexp(Min(kappa - nu,n - q),1) - psiexp(1 + Min(kappa - nu,n - q),1) + n*psiexp(1 + Min(kappa - nu,n - q),1) - q*psiexp(1 + Min(kappa - nu,n - q),1)) + vEG2G2G2cAD*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(j - kappa,1)*(-psiexp(1 + Max(0,-j + kappa),1) + n*psiexp(1 + Max(0,-j + kappa),1) - q*psiexp(1 + Max(0,-j + kappa),1) + Max(0,-j + kappa)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,-j + kappa),1) - n*psiexp(2 + Max(0,-j + kappa),1) + q*psiexp(2 + Max(0,-j + kappa),1) - n*psiexp(Min(kappa - nu,n - q),1) + q*psiexp(Min(kappa - nu,n - q),1) - Min(kappa - nu,n - q)*(-1 + psiexp(1,1))*psiexp(Min(kappa - nu,n - q),1) + psiexp(1 + Min(kappa - nu,n - q),1) + n*psiexp(1 + Min(kappa - nu,n - q),1) - q*psiexp(1 + Min(kappa - nu,n - q),1)) + vEG2G2G2cAWBD*b[kappa]*Power2(1/(psiexp(1,1) - psiexp(1,2)))*psiexp(j - kappa,2)*psiexp(kappa - nu - Max(0,-j + kappa) - Min(kappa - nu,n - q),1)*(Max(0,-j + kappa)*(-psiexp(1,1) + psiexp(1,2))*psiexp(1 + Max(0,-j + kappa),2)*psiexp(Min(kappa - nu,n - q),1) - n*psiexp(2 + Max(0,-j + kappa),2)*psiexp(Min(kappa - nu,n - q),1) + q*psiexp(2 + Max(0,-j + kappa),2)*psiexp(Min(kappa - nu,n - q),1) + Min(kappa - nu,n - q)*(psiexp(1,1) - psiexp(1,2))*psiexp(1 + Max(0,-j + kappa),1)*psiexp(Min(kappa - nu,n - q),2) - n*psiexp(2 + Max(0,-j + kappa),1)*psiexp(Min(kappa - nu,n - q),2) + q*psiexp(2 + Max(0,-j + kappa),1)*psiexp(Min(kappa - nu,n - q),2) - psiexp(1 + Max(0,-j + kappa),2)*psiexp(1 + Min(kappa - nu,n - q),1) + n*psiexp(1 + Max(0,-j + kappa),2)*psiexp(1 + Min(kappa - nu,n - q),1) - q*psiexp(1 + Max(0,-j + kappa),2)*psiexp(1 + Min(kappa - nu,n - q),1) + psiexp(1 + Max(0,-j + kappa),1)*psiexp(1 + Min(kappa - nu,n - q),2) + n*psiexp(1 + Max(0,-j + kappa),1)*psiexp(1 + Min(kappa - nu,n - q),2) - q*psiexp(1 + Max(0,-j + kappa),1)*psiexp(1 + Min(kappa - nu,n - q),2));
					}
					sumInuP1 = sumIkappa;
				}
				sumInu += sumInuP1*b[nu];
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += MMjkP1/(n - q);
	}
}

for (int j = 1; j <= q; j++)
{
	for (int k = 0; k <= 0; k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = 0; nu <= -1 + j; nu++)
			{
				scalar sumInuP1;
				{
					scalar sumIkappa = 0;
					for (int kappa = 0; kappa <= Min(-2 + n + nu - q,q); kappa++)
					{
						sumIkappa += (vEG2G2G2c*b[kappa]*(-n + q - 2*n*q + (1 - 2*n + 2*q)*Max(0,kappa - nu) + Power2(n) + Power2(q) + Power2(Max(0,kappa - nu))))/2. + (vEG2G2G2cAD*b[kappa]*(-n + q - 2*n*q + (1 - 2*n + 2*q)*Max(0,kappa - nu) + Power2(n) + Power2(q) + Power2(Max(0,kappa - nu)))*psiexp(j - nu,1))/2. + vEG2G2G2cAWAD*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(j - kappa,1)*(psiexp(1 + n - q,1) + (-1 + n - q)*psiexp(1 + Max(0,kappa - nu),1) + Max(0,kappa - nu)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,kappa - nu),1) + (-n + q)*psiexp(2 + Max(0,kappa - nu),1)) + vEG2G2G2cAW*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-kappa + nu,1)*(psiexp(1 + n - q,1) + (-1 + n - q)*psiexp(1 + Max(0,kappa - nu),1) + Max(0,kappa - nu)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,kappa - nu),1) + (-n + q)*psiexp(2 + Max(0,kappa - nu),1)) + vEG2G2G2cAWBD*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(j - nu,2)*psiexp(-kappa + nu,1)*(psiexp(1 + n - q,1) + (-1 + n - q)*psiexp(1 + Max(0,kappa - nu),1) + Max(0,kappa - nu)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,kappa - nu),1) + (-n + q)*psiexp(2 + Max(0,kappa - nu),1));
					}
					sumInuP1 = sumIkappa;
				}
				sumInu += sumInuP1*b[nu];
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += MMjkP1/(n - q);
	}
}

for (int j = 1; j <= -1 + q; j++)
{
	for (int k = 0; k <= 0; k++)
	{
		scalar MMjkP1;
		{
			scalar sumIkappa = 0;
			for (int kappa = Max(0,1 + j); kappa <= Min(-1 + j + n - q,q); kappa++)
			{
				sumIkappa += (j - kappa + n - q)*b[kappa];
			}
			MMjkP1 = sumIkappa;
		}
		scalar MMjkP2;
		{
			scalar sumInu = 0;
			for (int nu = 0; nu <= -1 + j; nu++)
			{
				sumInu += vEG2G4c*b[nu] + vEG2G4cAW*b[nu]*psiexp(j - nu,1);
			}
			MMjkP2 = sumInu;
		}
		MM(j, k) += (MMjkP1*MMjkP2)/(n - q);
	}
}

for (int j = 1; j <= q; j++)
{
	for (int k = 0; k <= 0; k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = 0; nu <= Min(-1 + j,-1 + q); nu++)
			{
				scalar sumInuP1;
				{
					scalar sumIkappa = 0;
					for (int kappa = Max(0,1 + nu); kappa <= Min(-1 + n + nu - q,q); kappa++)
					{
						sumIkappa += (-kappa + n + nu - q)*b[kappa];
					}
					sumInuP1 = sumIkappa;
				}
				sumInu += sumInuP1*b[nu]*(vEG4G2c + vEG4G2cAD*psiexp(j - nu,1) + vEG4G2cBD*psiexp(j - nu,2));
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += MMjkP1/(n - q);
	}
}

for (int j = 1; j <= -2 + q; j++)
{
	for (int k = 0; k <= 0; k++)
	{
		scalar MMjkP1;
		{
			scalar sumIkappa = 0;
			for (int kappa = Max(0,2 + j); kappa <= q; kappa++)
			{
				sumIkappa += (vEG4G2c*(2*n - 2*q - Min(-j + kappa,n - q))*(-1 + Min(-j + kappa,n - q))*b[kappa])/2. + vEG4G2cAD*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-j + kappa - Min(-j + kappa,n - q),1)*((1 + n - q)*psiexp(1,1) + Min(-j + kappa,n - q)*(-1 + psiexp(1,1))*psiexp(1,1) + (-n + q)*psiexp(2,1) + (-n + q)*psiexp(Min(-j + kappa,n - q),1) + (-1 + n - q)*psiexp(1 + Min(-j + kappa,n - q),1)) + vEG4G2cBD*b[kappa]*Power2(1/(-1 + psiexp(1,2)))*psiexp(-j + kappa - Min(-j + kappa,n - q),2)*((1 + n - q)*psiexp(1,2) + Min(-j + kappa,n - q)*(-1 + psiexp(1,2))*psiexp(1,2) + (-n + q)*psiexp(2,2) + (-n + q)*psiexp(Min(-j + kappa,n - q),2) + (-1 + n - q)*psiexp(1 + Min(-j + kappa,n - q),2));
			}
			MMjkP1 = sumIkappa;
		}
		MM(j, k) += (MMjkP1*b[j])/(n - q);
	}
}

for (int j = 1; j <= q; j++)
{
	for (int k = 0; k <= 0; k++)
	{
		scalar MMjkP1;
		{
			scalar sumIkappa = 0;
			for (int kappa = 0; kappa <= Min(-2 + j + n - q,q); kappa++)
			{
				sumIkappa += (vEG2G4c*b[kappa]*(-n + q - 2*n*q + (1 - 2*n + 2*q)*Max(0,-j + kappa) + Power2(n) + Power2(q) + Power2(Max(0,-j + kappa))))/2. + vEG2G4cAW*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(j - kappa,1)*(psiexp(1 + n - q,1) + (-1 + n - q)*psiexp(1 + Max(0,-j + kappa),1) + Max(0,-j + kappa)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,-j + kappa),1) + (-n + q)*psiexp(2 + Max(0,-j + kappa),1));
			}
			MMjkP1 = sumIkappa;
		}
		MM(j, k) += (MMjkP1*b[j])/(n - q);
	}
}

for (int j = 1; j <= -1 + q; j++)
{
	for (int k = 0; k <= 0; k++)
	{
		scalar MMjkP1;
		{
			scalar sumIkappa = 0;
			for (int kappa = Max(0,1 + j); kappa <= Min(-1 + j + n - q,q); kappa++)
			{
				sumIkappa += (j - kappa + n - q)*b[kappa];
			}
			MMjkP1 = sumIkappa;
		}
		MM(j, k) += (MMjkP1*vEG6c*b[j])/(n - q);
	}
}

for (int j = 1; j <= -2 + q; j++)
{
	for (int k = 0; k <= 0; k++)
	{
		scalar MMjkP1;
		{
			scalar sumIkappa = 0;
			for (int kappa = Max(0,2 + j); kappa <= q; kappa++)
			{
				sumIkappa += (vEG2G2c*(2*n - 2*q - Min(-j + kappa,n - q))*(-1 + Min(-j + kappa,n - q))*b[kappa])/2. + vEG2G2cAW*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-j + kappa - Min(-j + kappa,n - q),1)*((1 + n - q)*psiexp(1,1) + Min(-j + kappa,n - q)*(-1 + psiexp(1,1))*psiexp(1,1) + (-n + q)*psiexp(2,1) + (-n + q)*psiexp(Min(-j + kappa,n - q),1) + (-1 + n - q)*psiexp(1 + Min(-j + kappa,n - q),1));
			}
			MMjkP1 = sumIkappa;
		}
		MM(j, k) += (azero*MMjkP1)/(n - q);
	}
}

for (int j = 1; j <= q; j++)
{
	for (int k = 0; k <= 0; k++)
	{
		scalar MMjkP1;
		{
			scalar sumIkappa = 0;
			for (int kappa = 0; kappa <= Min(-2 + j + n - q,q); kappa++)
			{
				sumIkappa += (vEG2G2c*b[kappa]*(-n + q - 2*n*q + (1 - 2*n + 2*q)*Max(0,-j + kappa) + Power2(n) + Power2(q) + Power2(Max(0,-j + kappa))))/2. + vEG2G2cAW*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(j - kappa,1)*(psiexp(1 + n - q,1) + (-1 + n - q)*psiexp(1 + Max(0,-j + kappa),1) + Max(0,-j + kappa)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,-j + kappa),1) + (-n + q)*psiexp(2 + Max(0,-j + kappa),1));
			}
			MMjkP1 = sumIkappa;
		}
		MM(j, k) += (azero*MMjkP1)/(n - q);
	}
}

for (int j = 1; j <= -1 + q; j++)
{
	for (int k = 0; k <= 0; k++)
	{
		scalar MMjkP1;
		{
			scalar sumIkappa = 0;
			for (int kappa = Max(0,1 + j); kappa <= Min(-1 + j + n - q,q); kappa++)
			{
				sumIkappa += (j - kappa + n - q)*b[kappa];
			}
			MMjkP1 = sumIkappa;
		}
		MM(j, k) += (azero*MMjkP1*vEG4c)/(n - q);
	}
}

for (int j = 1; j <= -1 + q; j++)
{
	for (int k = 0; k <= 0; k++)
	{
		scalar MMjkP1;
		{
			scalar sumIkappa = 0;
			for (int kappa = 1 + j; kappa <= q; kappa++)
			{
				sumIkappa += vEG2G2c*b[kappa] + vEG2G2cAW*b[kappa]*psiexp(-j + kappa,1);
			}
			MMjkP1 = sumIkappa;
		}
		MM(j, k) += (azero*MMjkP1*(-n + q - 2*n*q + Power2(n) + Power2(q)))/(2.*(n - q));
	}
}

for (int j = 1; j <= q; j++)
{
	for (int k = 0; k <= 0; k++)
	{
		scalar MMjkP1;
		{
			scalar sumIkappa = 0;
			for (int kappa = 0; kappa <= -1 + j; kappa++)
			{
				sumIkappa += vEG2G2c*b[kappa] + vEG2G2cAW*b[kappa]*psiexp(j - kappa,1);
			}
			MMjkP1 = sumIkappa;
		}
		MM(j, k) += (azero*MMjkP1*(-n + q - 2*n*q + Power2(n) + Power2(q)))/(2.*(n - q));
	}
}

for (int j = 1; j <= q; j++)
{
	for (int k = 0; k <= 0; k++)
	{
		MM(j, k) += (azero*(-1 + n - q)*vEG4c*b[j])/2.;
	}
}

for (int j = 1; j <= q; j++)
{
	for (int k = 0; k <= 0; k++)
	{
		MM(j, k) += ((-1 + n - q)*vEG2c*Power2(azero))/2.;
	}
}

for (int j = 0; j <= 0; j++)
{
	for (int k = 2; k <= -1 + q; k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = 0; nu <= Min(-2 + k,q); nu++)
			{
				scalar sumInuP1;
				{
					scalar sumIkappa = 0;
					for (int kappa = 1 + k; kappa <= q; kappa++)
					{
						sumIkappa += (vEG2G2G2c*(2*n - 2*q - Min(k - nu,n - q))*(-1 + Min(k - nu,n - q))*b[kappa])/2. + (vEG2G2G2cAD*(2*n - 2*q - Min(k - nu,n - q))*(-1 + Min(k - nu,n - q))*b[kappa]*psiexp(-k + kappa,1))/2. + vEG2G2G2cAW*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(k - nu - Min(k - nu,n - q),1)*((1 + n - q)*psiexp(1,1) + Min(k - nu,n - q)*(-1 + psiexp(1,1))*psiexp(1,1) + (-n + q)*psiexp(2,1) + (-n + q)*psiexp(Min(k - nu,n - q),1) + (-1 + n - q)*psiexp(1 + Min(k - nu,n - q),1)) + vEG2G2G2cAWBD*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-k + kappa,2)*psiexp(k - nu - Min(k - nu,n - q),1)*((1 + n - q)*psiexp(1,1) + Min(k - nu,n - q)*(-1 + psiexp(1,1))*psiexp(1,1) + (-n + q)*psiexp(2,1) + (-n + q)*psiexp(Min(k - nu,n - q),1) + (-1 + n - q)*psiexp(1 + Min(k - nu,n - q),1)) + vEG2G2G2cAWAD*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(kappa - nu - Min(k - nu,n - q),1)*((1 + n - q)*psiexp(1,1) + Min(k - nu,n - q)*(-1 + psiexp(1,1))*psiexp(1,1) + (-n + q)*psiexp(2,1) + (-n + q)*psiexp(Min(k - nu,n - q),1) + (-1 + n - q)*psiexp(1 + Min(k - nu,n - q),1));
					}
					sumInuP1 = sumIkappa;
				}
				sumInu += sumInuP1*b[nu];
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += MMjkP1/(n - q);
	}
}

for (int j = 0; j <= 0; j++)
{
	for (int k = 3; k <= q; k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = 0; nu <= Min(-3 + k,q); nu++)
			{
				scalar sumInuP1;
				{
					scalar sumIkappa = 0;
					for (int kappa = Max(0,2 + nu); kappa <= -1 + k; kappa++)
					{
						sumIkappa += (vEG2G2G2c*(2*n - 2*q - Min(kappa - nu,n - q))*(-1 + Min(kappa - nu,n - q))*b[kappa])/2. + (vEG2G2G2cAD*(2*n - 2*q - Min(kappa - nu,n - q))*(-1 + Min(kappa - nu,n - q))*b[kappa]*psiexp(k - kappa,1))/2. + vEG2G2G2cAWAD*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(k - nu - Min(kappa - nu,n - q),1)*((1 + n - q)*psiexp(1,1) + Min(kappa - nu,n - q)*(-1 + psiexp(1,1))*psiexp(1,1) + (-n + q)*psiexp(2,1) + (-n + q)*psiexp(Min(kappa - nu,n - q),1) + (-1 + n - q)*psiexp(1 + Min(kappa - nu,n - q),1)) + vEG2G2G2cAW*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(kappa - nu - Min(kappa - nu,n - q),1)*((1 + n - q)*psiexp(1,1) + Min(kappa - nu,n - q)*(-1 + psiexp(1,1))*psiexp(1,1) + (-n + q)*psiexp(2,1) + (-n + q)*psiexp(Min(kappa - nu,n - q),1) + (-1 + n - q)*psiexp(1 + Min(kappa - nu,n - q),1)) + vEG2G2G2cAWBD*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(k - kappa,2)*psiexp(kappa - nu - Min(kappa - nu,n - q),1)*((1 + n - q)*psiexp(1,1) + Min(kappa - nu,n - q)*(-1 + psiexp(1,1))*psiexp(1,1) + (-n + q)*psiexp(2,1) + (-n + q)*psiexp(Min(kappa - nu,n - q),1) + (-1 + n - q)*psiexp(1 + Min(kappa - nu,n - q),1));
					}
					sumInuP1 = sumIkappa;
				}
				sumInu += sumInuP1*b[nu];
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += MMjkP1/(n - q);
	}
}

for (int j = 0; j <= 0; j++)
{
	for (int k = 2; k <= q; k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = 0; nu <= Min(-2 + k,q); nu++)
			{
				scalar sumInuP1;
				{
					scalar sumIkappa = 0;
					for (int kappa = 0; kappa <= Min(-2 + k,-2 + n + nu - q); kappa++)
					{
						sumIkappa += -(vEG2G2G2c*(2*n - 2*q - Max(0,kappa - nu) - Min(k - nu,n - q))*(1 + Max(0,kappa - nu) - Min(k - nu,n - q))*b[kappa])/2. - (vEG2G2G2cAWAD*(2*n - 2*q - Max(0,kappa - nu) - Min(k - nu,n - q))*(1 + Max(0,kappa - nu) - Min(k - nu,n - q))*b[kappa]*psiexp(k - kappa,1))/2. + vEG2G2G2cAD*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(k - nu - Max(0,kappa - nu) - Min(k - nu,n - q),1)*(psiexp(1 + Max(0,kappa - nu),1) + n*psiexp(1 + Max(0,kappa - nu),1) - q*psiexp(1 + Max(0,kappa - nu),1) + Min(k - nu,n - q)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,kappa - nu),1) - n*psiexp(2 + Max(0,kappa - nu),1) + q*psiexp(2 + Max(0,kappa - nu),1) - n*psiexp(Min(k - nu,n - q),1) + q*psiexp(Min(k - nu,n - q),1) - Max(0,kappa - nu)*(-1 + psiexp(1,1))*psiexp(Min(k - nu,n - q),1) - psiexp(1 + Min(k - nu,n - q),1) + n*psiexp(1 + Min(k - nu,n - q),1) - q*psiexp(1 + Min(k - nu,n - q),1)) + vEG2G2G2cAW*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-kappa + nu,1)*(-psiexp(1 + Max(0,kappa - nu),1) + n*psiexp(1 + Max(0,kappa - nu),1) - q*psiexp(1 + Max(0,kappa - nu),1) + Max(0,kappa - nu)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,kappa - nu),1) - n*psiexp(2 + Max(0,kappa - nu),1) + q*psiexp(2 + Max(0,kappa - nu),1) - n*psiexp(Min(k - nu,n - q),1) + q*psiexp(Min(k - nu,n - q),1) - Min(k - nu,n - q)*(-1 + psiexp(1,1))*psiexp(Min(k - nu,n - q),1) + psiexp(1 + Min(k - nu,n - q),1) + n*psiexp(1 + Min(k - nu,n - q),1) - q*psiexp(1 + Min(k - nu,n - q),1)) + vEG2G2G2cAWBD*b[kappa]*Power2(1/(psiexp(1,1) - psiexp(1,2)))*psiexp(-kappa + nu,1)*psiexp(k - nu - Max(0,kappa - nu) - Min(k - nu,n - q),2)*(Min(k - nu,n - q)*(-psiexp(1,1) + psiexp(1,2))*psiexp(1 + Max(0,kappa - nu),2)*psiexp(Min(k - nu,n - q),1) - n*psiexp(2 + Max(0,kappa - nu),2)*psiexp(Min(k - nu,n - q),1) + q*psiexp(2 + Max(0,kappa - nu),2)*psiexp(Min(k - nu,n - q),1) + Max(0,kappa - nu)*(psiexp(1,1) - psiexp(1,2))*psiexp(1 + Max(0,kappa - nu),1)*psiexp(Min(k - nu,n - q),2) - n*psiexp(2 + Max(0,kappa - nu),1)*psiexp(Min(k - nu,n - q),2) + q*psiexp(2 + Max(0,kappa - nu),1)*psiexp(Min(k - nu,n - q),2) + psiexp(1 + Max(0,kappa - nu),2)*psiexp(1 + Min(k - nu,n - q),1) + n*psiexp(1 + Max(0,kappa - nu),2)*psiexp(1 + Min(k - nu,n - q),1) - q*psiexp(1 + Max(0,kappa - nu),2)*psiexp(1 + Min(k - nu,n - q),1) - psiexp(1 + Max(0,kappa - nu),1)*psiexp(1 + Min(k - nu,n - q),2) + n*psiexp(1 + Max(0,kappa - nu),1)*psiexp(1 + Min(k - nu,n - q),2) - q*psiexp(1 + Max(0,kappa - nu),1)*psiexp(1 + Min(k - nu,n - q),2));
					}
					sumInuP1 = sumIkappa;
				}
				sumInu += sumInuP1*b[nu];
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += MMjkP1/(n - q);
	}
}

for (int j = 0; j <= 0; j++)
{
	for (int k = 2; k <= q; k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = 0; nu <= Min(-2 + k,q); nu++)
			{
				sumInu += (vEG2G4c*(2*n - 2*q - Min(k - nu,n - q))*(-1 + Min(k - nu,n - q))*b[nu])/2. + vEG2G4cAW*b[nu]*Power2(1/(-1 + psiexp(1,1)))*psiexp(k - nu - Min(k - nu,n - q),1)*((1 + n - q)*psiexp(1,1) + Min(k - nu,n - q)*(-1 + psiexp(1,1))*psiexp(1,1) + (-n + q)*psiexp(2,1) + (-n + q)*psiexp(Min(k - nu,n - q),1) + (-1 + n - q)*psiexp(1 + Min(k - nu,n - q),1));
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += (MMjkP1*b[k])/(n - q);
	}
}

for (int j = 0; j <= 0; j++)
{
	for (int k = 2; k <= q; k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = 0; nu <= Min(-2 + k,q); nu++)
			{
				scalar sumInuP1;
				{
					scalar sumIkappa = 0;
					for (int kappa = Max(0,1 + nu); kappa <= Min(-1 + k,-1 + n + nu - q); kappa++)
					{
						sumIkappa += (-kappa + n + nu - q)*vEG4G2c*b[kappa] - (kappa - n - nu + q)*vEG4G2cAD*b[kappa]*psiexp(k - kappa,1) - (kappa - n - nu + q)*vEG4G2cBD*b[kappa]*psiexp(k - kappa,2);
					}
					sumInuP1 = sumIkappa;
				}
				sumInu += sumInuP1*b[nu];
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += MMjkP1/(n - q);
	}
}

for (int j = 0; j <= 0; j++)
{
	for (int k = 1; k <= -2 + q; k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = Max(0,2 + k - n + q); nu <= -2 + q; nu++)
			{
				scalar sumInuP1;
				{
					scalar sumIkappa = 0;
					for (int kappa = Max(2 + k,2 + nu); kappa <= q; kappa++)
					{
						sumIkappa += -(vEG2G2G2c*(2*n - 2*q - Max(0,k - nu) - Min(kappa - nu,n - q))*(1 + Max(0,k - nu) - Min(kappa - nu,n - q))*b[kappa])/2. - (vEG2G2G2cAWAD*(2*n - 2*q - Max(0,k - nu) - Min(kappa - nu,n - q))*(1 + Max(0,k - nu) - Min(kappa - nu,n - q))*b[kappa]*psiexp(-k + kappa,1))/2. + vEG2G2G2cAD*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(kappa - nu - Max(0,k - nu) - Min(kappa - nu,n - q),1)*(psiexp(1 + Max(0,k - nu),1) + n*psiexp(1 + Max(0,k - nu),1) - q*psiexp(1 + Max(0,k - nu),1) + Min(kappa - nu,n - q)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,k - nu),1) - n*psiexp(2 + Max(0,k - nu),1) + q*psiexp(2 + Max(0,k - nu),1) - n*psiexp(Min(kappa - nu,n - q),1) + q*psiexp(Min(kappa - nu,n - q),1) - Max(0,k - nu)*(-1 + psiexp(1,1))*psiexp(Min(kappa - nu,n - q),1) - psiexp(1 + Min(kappa - nu,n - q),1) + n*psiexp(1 + Min(kappa - nu,n - q),1) - q*psiexp(1 + Min(kappa - nu,n - q),1)) + vEG2G2G2cAW*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-k + nu,1)*(-psiexp(1 + Max(0,k - nu),1) + n*psiexp(1 + Max(0,k - nu),1) - q*psiexp(1 + Max(0,k - nu),1) + Max(0,k - nu)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,k - nu),1) - n*psiexp(2 + Max(0,k - nu),1) + q*psiexp(2 + Max(0,k - nu),1) - n*psiexp(Min(kappa - nu,n - q),1) + q*psiexp(Min(kappa - nu,n - q),1) - Min(kappa - nu,n - q)*(-1 + psiexp(1,1))*psiexp(Min(kappa - nu,n - q),1) + psiexp(1 + Min(kappa - nu,n - q),1) + n*psiexp(1 + Min(kappa - nu,n - q),1) - q*psiexp(1 + Min(kappa - nu,n - q),1)) + vEG2G2G2cAWBD*b[kappa]*Power2(1/(psiexp(1,1) - psiexp(1,2)))*psiexp(-k + nu,1)*psiexp(kappa - nu - Max(0,k - nu) - Min(kappa - nu,n - q),2)*(Min(kappa - nu,n - q)*(-psiexp(1,1) + psiexp(1,2))*psiexp(1 + Max(0,k - nu),2)*psiexp(Min(kappa - nu,n - q),1) - n*psiexp(2 + Max(0,k - nu),2)*psiexp(Min(kappa - nu,n - q),1) + q*psiexp(2 + Max(0,k - nu),2)*psiexp(Min(kappa - nu,n - q),1) + Max(0,k - nu)*(psiexp(1,1) - psiexp(1,2))*psiexp(1 + Max(0,k - nu),1)*psiexp(Min(kappa - nu,n - q),2) - n*psiexp(2 + Max(0,k - nu),1)*psiexp(Min(kappa - nu,n - q),2) + q*psiexp(2 + Max(0,k - nu),1)*psiexp(Min(kappa - nu,n - q),2) + psiexp(1 + Max(0,k - nu),2)*psiexp(1 + Min(kappa - nu,n - q),1) + n*psiexp(1 + Max(0,k - nu),2)*psiexp(1 + Min(kappa - nu,n - q),1) - q*psiexp(1 + Max(0,k - nu),2)*psiexp(1 + Min(kappa - nu,n - q),1) - psiexp(1 + Max(0,k - nu),1)*psiexp(1 + Min(kappa - nu,n - q),2) + n*psiexp(1 + Max(0,k - nu),1)*psiexp(1 + Min(kappa - nu,n - q),2) - q*psiexp(1 + Max(0,k - nu),1)*psiexp(1 + Min(kappa - nu,n - q),2));
					}
					sumInuP1 = sumIkappa;
				}
				sumInu += sumInuP1*b[nu];
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += MMjkP1/(n - q);
	}
}

for (int j = 0; j <= 0; j++)
{
	for (int k = 1; k <= -1 + q; k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = Max(0,3 + k - n + q); nu <= q; nu++)
			{
				scalar sumInuP1;
				{
					scalar sumIkappa = 0;
					for (int kappa = 1 + k; kappa <= Min(-2 + n + nu - q,q); kappa++)
					{
						sumIkappa += (vEG2G2G2c*b[kappa]*(-n + q - 2*n*q + (1 - 2*n + 2*q)*Max(0,kappa - nu) + Power2(n) + Power2(q) + Power2(Max(0,kappa - nu))))/2. + (vEG2G2G2cAW*b[kappa]*(-n + q - 2*n*q + (1 - 2*n + 2*q)*Max(0,kappa - nu) + Power2(n) + Power2(q) + Power2(Max(0,kappa - nu)))*psiexp(-k + kappa,1))/2. + vEG2G2G2cAWAD*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-k + nu,1)*(psiexp(1 + n - q,1) + (-1 + n - q)*psiexp(1 + Max(0,kappa - nu),1) + Max(0,kappa - nu)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,kappa - nu),1) + (-n + q)*psiexp(2 + Max(0,kappa - nu),1)) + vEG2G2G2cAD*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-kappa + nu,1)*(psiexp(1 + n - q,1) + (-1 + n - q)*psiexp(1 + Max(0,kappa - nu),1) + Max(0,kappa - nu)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,kappa - nu),1) + (-n + q)*psiexp(2 + Max(0,kappa - nu),1)) + vEG2G2G2cAWBD*b[kappa]*Power2(1/(-1 + psiexp(1,2)))*psiexp(-k + kappa,1)*psiexp(-kappa + nu,2)*(psiexp(1 + n - q,2) + (-1 + n - q)*psiexp(1 + Max(0,kappa - nu),2) + Max(0,kappa - nu)*(-1 + psiexp(1,2))*psiexp(1 + Max(0,kappa - nu),2) + (-n + q)*psiexp(2 + Max(0,kappa - nu),2));
					}
					sumInuP1 = sumIkappa;
				}
				sumInu += sumInuP1*b[nu];
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += MMjkP1/(n - q);
	}
}

for (int j = 0; j <= 0; j++)
{
	for (int k = 1; k <= q; k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = Max(0,2 + k - n + q); nu <= q; nu++)
			{
				scalar sumInuP1;
				{
					scalar sumIkappa = 0;
					for (int kappa = 0; kappa <= -1 + k; kappa++)
					{
						sumIkappa += (vEG2G2G2c*b[kappa]*(-n + q - 2*n*q + (1 - 2*n + 2*q)*Max(0,k - nu) + Power2(n) + Power2(q) + Power2(Max(0,k - nu))))/2. + (vEG2G2G2cAW*b[kappa]*(-n + q - 2*n*q + (1 - 2*n + 2*q)*Max(0,k - nu) + Power2(n) + Power2(q) + Power2(Max(0,k - nu)))*psiexp(k - kappa,1))/2. + vEG2G2G2cAD*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-k + nu,1)*(psiexp(1 + n - q,1) + (-1 + n - q)*psiexp(1 + Max(0,k - nu),1) + Max(0,k - nu)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,k - nu),1) + (-n + q)*psiexp(2 + Max(0,k - nu),1)) + vEG2G2G2cAWAD*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-kappa + nu,1)*(psiexp(1 + n - q,1) + (-1 + n - q)*psiexp(1 + Max(0,k - nu),1) + Max(0,k - nu)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,k - nu),1) + (-n + q)*psiexp(2 + Max(0,k - nu),1)) + vEG2G2G2cAWBD*b[kappa]*Power2(1/(-1 + psiexp(1,2)))*psiexp(k - kappa,1)*psiexp(-k + nu,2)*(psiexp(1 + n - q,2) + (-1 + n - q)*psiexp(1 + Max(0,k - nu),2) + Max(0,k - nu)*(-1 + psiexp(1,2))*psiexp(1 + Max(0,k - nu),2) + (-n + q)*psiexp(2 + Max(0,k - nu),2));
					}
					sumInuP1 = sumIkappa;
				}
				sumInu += sumInuP1*b[nu];
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += MMjkP1/(n - q);
	}
}

for (int j = 0; j <= 0; j++)
{
	for (int k = 1; k <= -1 + q; k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = Max(0,2 + k - n + q); nu <= -1 + q; nu++)
			{
				scalar sumInuP1;
				{
					scalar sumIkappa = 0;
					for (int kappa = Max(1 + k,1 + nu); kappa <= Min(-1 + n + nu - q,q); kappa++)
					{
						sumIkappa += (-kappa + n + nu - q)*vEG2G4c*b[kappa] + (-kappa + n + nu - q)*vEG2G4cAW*b[kappa]*psiexp(-k + kappa,1);
					}
					sumInuP1 = sumIkappa;
				}
				sumInu += sumInuP1*b[nu];
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += MMjkP1/(n - q);
	}
}

for (int j = 0; j <= 0; j++)
{
	for (int k = 1; k <= q; k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = Max(0,2 + k - n + q); nu <= q; nu++)
			{
				sumInu += (vEG4G2c*b[nu]*(-n + q - 2*n*q + (1 - 2*n + 2*q)*Max(0,k - nu) + Power2(n) + Power2(q) + Power2(Max(0,k - nu))))/2. + vEG4G2cAD*b[nu]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-k + nu,1)*(psiexp(1 + n - q,1) + (-1 + n - q)*psiexp(1 + Max(0,k - nu),1) + Max(0,k - nu)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,k - nu),1) + (-n + q)*psiexp(2 + Max(0,k - nu),1)) + vEG4G2cBD*b[nu]*Power2(1/(-1 + psiexp(1,2)))*psiexp(-k + nu,2)*(psiexp(1 + n - q,2) + (-1 + n - q)*psiexp(1 + Max(0,k - nu),2) + Max(0,k - nu)*(-1 + psiexp(1,2))*psiexp(1 + Max(0,k - nu),2) + (-n + q)*psiexp(2 + Max(0,k - nu),2));
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += (MMjkP1*b[k])/(n - q);
	}
}

for (int j = 0; j <= 0; j++)
{
	for (int k = 1; k <= -1 + q; k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = Max(0,1 + k - n + q); nu <= Min(-1 + k,q); nu++)
			{
				scalar sumInuP1;
				{
					scalar sumIkappa = 0;
					for (int kappa = 1 + k; kappa <= q; kappa++)
					{
						sumIkappa += (-k + n + nu - q)*vEG4G2c*b[kappa] - (k - n - nu + q)*vEG4G2cAD*b[kappa]*psiexp(-k + kappa,1) - (k - n - nu + q)*vEG4G2cBD*b[kappa]*psiexp(-k + kappa,2);
					}
					sumInuP1 = sumIkappa;
				}
				sumInu += sumInuP1*b[nu];
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += MMjkP1/(n - q);
	}
}

for (int j = 0; j <= 0; j++)
{
	for (int k = 1; k <= q; k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = Max(0,1 + k - n + q); nu <= Min(-1 + k,q); nu++)
			{
				scalar sumInuP1;
				{
					scalar sumIkappa = 0;
					for (int kappa = 0; kappa <= -1 + k; kappa++)
					{
						sumIkappa += (-k + n + nu - q)*vEG2G4c*b[kappa] - (k - n - nu + q)*vEG2G4cAW*b[kappa]*psiexp(k - kappa,1);
					}
					sumInuP1 = sumIkappa;
				}
				sumInu += sumInuP1*b[nu];
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += MMjkP1/(n - q);
	}
}

for (int j = 0; j <= 0; j++)
{
	for (int k = 1; k <= q; k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = Max(0,1 + k - n + q); nu <= Min(-1 + k,q); nu++)
			{
				sumInu += (-k + n + nu - q)*b[nu];
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += (MMjkP1*vEG6c*b[k])/(n - q);
	}
}

for (int j = 0; j <= 0; j++)
{
	for (int k = 2; k <= q; k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = 0; nu <= Min(-2 + k,q); nu++)
			{
				sumInu += (vEG2G2c*(2*n - 2*q - Min(k - nu,n - q))*(-1 + Min(k - nu,n - q))*b[nu])/2. + vEG2G2cAW*b[nu]*Power2(1/(-1 + psiexp(1,1)))*psiexp(k - nu - Min(k - nu,n - q),1)*((1 + n - q)*psiexp(1,1) + Min(k - nu,n - q)*(-1 + psiexp(1,1))*psiexp(1,1) + (-n + q)*psiexp(2,1) + (-n + q)*psiexp(Min(k - nu,n - q),1) + (-1 + n - q)*psiexp(1 + Min(k - nu,n - q),1));
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += (azero*MMjkP1)/(n - q);
	}
}

for (int j = 0; j <= 0; j++)
{
	for (int k = 1; k <= q; k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = Max(0,2 + k - n + q); nu <= q; nu++)
			{
				sumInu += (vEG2G2c*b[nu]*(-n + q - 2*n*q + (1 - 2*n + 2*q)*Max(0,k - nu) + Power2(n) + Power2(q) + Power2(Max(0,k - nu))))/2. + vEG2G2cAW*b[nu]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-k + nu,1)*(psiexp(1 + n - q,1) + (-1 + n - q)*psiexp(1 + Max(0,k - nu),1) + Max(0,k - nu)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,k - nu),1) + (-n + q)*psiexp(2 + Max(0,k - nu),1));
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += (azero*MMjkP1)/(n - q);
	}
}

for (int j = 0; j <= 0; j++)
{
	for (int k = 1; k <= q; k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = Max(0,1 + k - n + q); nu <= Min(-1 + k,q); nu++)
			{
				sumInu += (-k + n + nu - q)*b[nu];
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += (azero*MMjkP1*vEG4c)/(n - q);
	}
}

for (int j = 0; j <= 0; j++)
{
	for (int k = 1; k <= -1 + q; k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = 1 + k; nu <= q; nu++)
			{
				sumInu += vEG2G2c*b[nu] + vEG2G2cAW*b[nu]*psiexp(-k + nu,1);
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += (azero*MMjkP1*(-n + q - 2*n*q + Power2(n) + Power2(q)))/(2.*(n - q));
	}
}

for (int j = 0; j <= 0; j++)
{
	for (int k = 1; k <= q; k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = 0; nu <= -1 + k; nu++)
			{
				sumInu += vEG2G2c*b[nu] + vEG2G2cAW*b[nu]*psiexp(k - nu,1);
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += (azero*MMjkP1*(-n + q - 2*n*q + Power2(n) + Power2(q)))/(2.*(n - q));
	}
}

for (int j = 0; j <= 0; j++)
{
	for (int k = 1; k <= q; k++)
	{
		MM(j, k) += (azero*(-1 + n - q)*vEG4c*b[k])/2.;
	}
}

for (int j = 0; j <= 0; j++)
{
	for (int k = 1; k <= q; k++)
	{
		MM(j, k) += ((-1 + n - q)*vEG2c*Power2(azero))/2.;
	}
}

for (int j = 1; j <= -2 + q; j++)
{
	for (int k = Max(1,2 + j); k <= q; k++)
	{
		MM(j, k) += -((1 + j - k)*(j - k + 2*n - 2*q)*vEG2G2c*Power2(azero))/(2.*(n - q)) + (vEG2G2cAW*Power2(azero)*Power2(1/(-1 + psiexp(1,1)))*(-n + q + (-1 + n - q)*psiexp(1,1) + (1 + j - k + n - q)*psiexp(1 + j - k,1) + (-j + k - n + q)*psiexp(2 + j - k,1))*psiexp(-j + k,1))/(n - q);
	}
}

for (int j = 1; j <= q; j++)
{
	for (int k = 1; k <= Min(-2 + j + n - q,q); k++)
	{
		MM(j, k) += (vEG2G2c*Power2(azero)*(-n + q - 2*n*q + (1 - 2*n + 2*q)*Max(0,-j + k) + Power2(n) + Power2(q) + Power2(Max(0,-j + k))))/(2.*(n - q)) + (vEG2G2cAW*Power2(azero)*Power2(1/(-1 + psiexp(1,1)))*psiexp(j - k,1)*(psiexp(1 + n - q,1) + (-1 + n - q)*psiexp(1 + Max(0,-j + k),1) + Max(0,-j + k)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,-j + k),1) + (-n + q)*psiexp(2 + Max(0,-j + k),1)))/(n - q);
	}
}

for (int j = 1; j <= -1 + q; j++)
{
	for (int k = Max(1,1 + j); k <= q; k++)
	{
		MM(j, k) += ((j - k + n - q)*vEG4c*Power2(azero))/(n - q);
	}
}

for (int j = 1; j <= -3 + q; j++)
{
	for (int k = Max(1,2 + j); k <= -1 + q; k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = 1 + k; nu <= q; nu++)
			{
				sumInu += -((1 + j - k)*(j - k + 2*n - 2*q)*vEG2G2G2c*b[nu])/2. + vEG2G2G2cAW*b[nu]*Power2(1/(-1 + psiexp(1,1)))*(-n + q + (-1 + n - q)*psiexp(1,1) + (1 + j - k + n - q)*psiexp(1 + j - k,1) + (-j + k - n + q)*psiexp(2 + j - k,1))*psiexp(-j + k,1) + vEG2G2G2cAWAD*b[nu]*Power2(1/(-1 + psiexp(1,1)))*(-n + q + (-1 + n - q)*psiexp(1,1) + (1 + j - k + n - q)*psiexp(1 + j - k,1) + (-j + k - n + q)*psiexp(2 + j - k,1))*psiexp(-j + nu,1) - ((1 + j - k)*(j - k + 2*n - 2*q)*vEG2G2G2cAD*b[nu]*psiexp(-k + nu,1))/2. + vEG2G2G2cAWBD*b[nu]*Power2(1/(-1 + psiexp(1,1)))*(-n + q + (-1 + n - q)*psiexp(1,1) + (1 + j - k + n - q)*psiexp(1 + j - k,1) + (-j + k - n + q)*psiexp(2 + j - k,1))*psiexp(-j + k,1)*psiexp(-k + nu,2);
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += (azero*MMjkP1)/(n - q);
	}
}

for (int j = 1; j <= -3 + q; j++)
{
	for (int k = Max(1,3 + j); k <= q; k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = Max(0,2 + j); nu <= -1 + k; nu++)
			{
				sumInu += -((1 + j - nu)*(j + 2*n - nu - 2*q)*vEG2G2G2c*b[nu])/2. + vEG2G2G2cAWAD*b[nu]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-j + k,1)*(-n + q + (-1 + n - q)*psiexp(1,1) + (1 + j + n - nu - q)*psiexp(1 + j - nu,1) + (-j - n + nu + q)*psiexp(2 + j - nu,1)) - ((-1 - j + nu)*(-j - 2*n + nu + 2*q)*vEG2G2G2cAD*b[nu]*psiexp(k - nu,1))/2. + vEG2G2G2cAW*b[nu]*Power2(1/(-1 + psiexp(1,1)))*(-n + q + (-1 + n - q)*psiexp(1,1) + (1 + j + n - nu - q)*psiexp(1 + j - nu,1) + (-j - n + nu + q)*psiexp(2 + j - nu,1))*psiexp(-j + nu,1) + vEG2G2G2cAWBD*b[nu]*Power2(1/(-1 + psiexp(1,1)))*(-n + q + (-1 + n - q)*psiexp(1,1) + (1 + j + n - nu - q)*psiexp(1 + j - nu,1) + (-j - n + nu + q)*psiexp(2 + j - nu,1))*psiexp(k - nu,2)*psiexp(-j + nu,1);
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += (azero*MMjkP1)/(n - q);
	}
}

for (int j = 1; j <= -2 + q; j++)
{
	for (int k = Max(2,2 + j); k <= q; k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = 0; nu <= -2 + k; nu++)
			{
				sumInu += -(vEG2G2G2c*(j - k + 2*n - 2*q - Max(0,-j + nu))*(1 + j - k + Max(0,-j + nu))*b[nu])/2. - (vEG2G2G2cAWAD*(-1 - j + k - Max(0,-j + nu))*(-j + k - 2*n + 2*q + Max(0,-j + nu))*b[nu]*psiexp(k - nu,1))/2. - vEG2G2G2cAD*b[nu]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-Max(0,-j + nu),1)*((n - q)*psiexp(-j + k,1) + Max(0,-j + nu)*(-1 + psiexp(1,1))*psiexp(-j + k,1) + (1 - n + q)*psiexp(1 - j + k,1) + (-1 - j + k - n + q)*psiexp(1 + Max(0,-j + nu),1) + (j - k + n - q)*psiexp(2 + Max(0,-j + nu),1)) + vEG2G2G2cAW*b[nu]*Power2(1/(-1 + psiexp(1,1)))*psiexp(j - nu,1)*((-j + k - n + q)*psiexp(-j + k,1) + (1 + j - k + n - q)*psiexp(1 - j + k,1) + (-1 + n - q)*psiexp(1 + Max(0,-j + nu),1) + Max(0,-j + nu)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,-j + nu),1) + (-n + q)*psiexp(2 + Max(0,-j + nu),1)) - vEG2G2G2cAWBD*b[nu]*Power2(1/(psiexp(1,1) - psiexp(1,2)))*psiexp(j - nu,1)*psiexp(-Max(0,-j + nu),2)*(-(Max(0,-j + nu)*(psiexp(1,1) - psiexp(1,2))*psiexp(-j + k,2)*psiexp(1 + Max(0,-j + nu),1)) + (1 - n + q)*psiexp(1 - j + k,2)*psiexp(1 + Max(0,-j + nu),1) + (-1 - j + k - n + q)*psiexp(1 - j + k,1)*psiexp(1 + Max(0,-j + nu),2) + (n - q)*psiexp(-j + k,2)*psiexp(2 + Max(0,-j + nu),1) + (j - k + n - q)*psiexp(-j + k,1)*psiexp(2 + Max(0,-j + nu),2));
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += (azero*MMjkP1)/(n - q);
	}
}

for (int j = 1; j <= -2 + q; j++)
{
	for (int k = Max(1,2 + j); k <= q; k++)
	{
		MM(j, k) += -(azero*(1 + j - k)*(j - k + 2*n - 2*q)*vEG2G4c*b[k])/(2.*(n - q)) + (azero*vEG2G4cAW*b[k]*Power2(1/(-1 + psiexp(1,1)))*(-n + q + (-1 + n - q)*psiexp(1,1) + (1 + j - k + n - q)*psiexp(1 + j - k,1) + (-j + k - n + q)*psiexp(2 + j - k,1))*psiexp(-j + k,1))/(n - q);
	}
}

for (int j = 1; j <= -2 + q; j++)
{
	for (int k = Max(1,2 + j); k <= q; k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = Max(0,1 + j); nu <= -1 + k; nu++)
			{
				sumInu += (j + n - nu - q)*vEG4G2c*b[nu] + (j + n - nu - q)*vEG4G2cAD*b[nu]*psiexp(k - nu,1) + (j + n - nu - q)*vEG4G2cBD*b[nu]*psiexp(k - nu,2);
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += (azero*MMjkP1)/(n - q);
	}
}

for (int j = 1; j <= -2 + q; j++)
{
	for (int k = 1; k <= -2 + q; k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = Max(2 + j,2 + k); nu <= q; nu++)
			{
				sumInu += -(vEG2G2G2c*(j + 2*n - nu - 2*q - Max(0,-j + k))*(1 + j - nu + Max(0,-j + k))*b[nu])/2. - (vEG2G2G2cAWAD*(j + 2*n - nu - 2*q - Max(0,-j + k))*(1 + j - nu + Max(0,-j + k))*b[nu]*psiexp(-k + nu,1))/2. - vEG2G2G2cAD*b[nu]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-Max(0,-j + k),1)*((n - q)*psiexp(-j + nu,1) + Max(0,-j + k)*(-1 + psiexp(1,1))*psiexp(-j + nu,1) + (1 - n + q)*psiexp(1 - j + nu,1) + (-1 - j - n + nu + q)*psiexp(1 + Max(0,-j + k),1) + (j + n - nu - q)*psiexp(2 + Max(0,-j + k),1)) + vEG2G2G2cAW*b[nu]*Power2(1/(-1 + psiexp(1,1)))*psiexp(j - k,1)*((-j - n + nu + q)*psiexp(-j + nu,1) + (1 + j + n - nu - q)*psiexp(1 - j + nu,1) + (-1 + n - q)*psiexp(1 + Max(0,-j + k),1) + Max(0,-j + k)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,-j + k),1) + (-n + q)*psiexp(2 + Max(0,-j + k),1)) - vEG2G2G2cAWBD*b[nu]*Power2(1/(psiexp(1,1) - psiexp(1,2)))*psiexp(j - k,1)*psiexp(-Max(0,-j + k),2)*(-(Max(0,-j + k)*(psiexp(1,1) - psiexp(1,2))*psiexp(-j + nu,2)*psiexp(1 + Max(0,-j + k),1)) + (1 - n + q)*psiexp(1 - j + nu,2)*psiexp(1 + Max(0,-j + k),1) + (-1 - j - n + nu + q)*psiexp(1 - j + nu,1)*psiexp(1 + Max(0,-j + k),2) + (n - q)*psiexp(-j + nu,2)*psiexp(2 + Max(0,-j + k),1) + (j + n - nu - q)*psiexp(-j + nu,1)*psiexp(2 + Max(0,-j + k),2));
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += (azero*MMjkP1)/(n - q);
	}
}

for (int j = 1; j <= q; j++)
{
	for (int k = 1; k <= Min(-3 + j + n - q,-1 + q); k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = 1 + k; nu <= Min(-2 + j + n - q,q); nu++)
			{
				sumInu += (vEG2G2G2c*b[nu]*(-n + q - 2*n*q + (1 - 2*n + 2*q)*Max(0,-j + nu) + Power2(n) + Power2(q) + Power2(Max(0,-j + nu))))/2. + (vEG2G2G2cAW*b[nu]*(-n + q - 2*n*q + (1 - 2*n + 2*q)*Max(0,-j + nu) + Power2(n) + Power2(q) + Power2(Max(0,-j + nu)))*psiexp(-k + nu,1))/2. + vEG2G2G2cAWAD*b[nu]*Power2(1/(-1 + psiexp(1,1)))*psiexp(j - k,1)*(psiexp(1 + n - q,1) + (-1 + n - q)*psiexp(1 + Max(0,-j + nu),1) + Max(0,-j + nu)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,-j + nu),1) + (-n + q)*psiexp(2 + Max(0,-j + nu),1)) + vEG2G2G2cAD*b[nu]*Power2(1/(-1 + psiexp(1,1)))*psiexp(j - nu,1)*(psiexp(1 + n - q,1) + (-1 + n - q)*psiexp(1 + Max(0,-j + nu),1) + Max(0,-j + nu)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,-j + nu),1) + (-n + q)*psiexp(2 + Max(0,-j + nu),1)) + vEG2G2G2cAWBD*b[nu]*Power2(1/(-1 + psiexp(1,2)))*psiexp(j - nu,2)*psiexp(-k + nu,1)*(psiexp(1 + n - q,2) + (-1 + n - q)*psiexp(1 + Max(0,-j + nu),2) + Max(0,-j + nu)*(-1 + psiexp(1,2))*psiexp(1 + Max(0,-j + nu),2) + (-n + q)*psiexp(2 + Max(0,-j + nu),2));
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += (azero*MMjkP1)/(n - q);
	}
}

for (int j = 1; j <= q; j++)
{
	for (int k = 1; k <= Min(-2 + j + n - q,q); k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = 0; nu <= -1 + k; nu++)
			{
				sumInu += (vEG2G2G2c*b[nu]*(-n + q - 2*n*q + (1 - 2*n + 2*q)*Max(0,-j + k) + Power2(n) + Power2(q) + Power2(Max(0,-j + k))))/2. + (vEG2G2G2cAW*b[nu]*(-n + q - 2*n*q + (1 - 2*n + 2*q)*Max(0,-j + k) + Power2(n) + Power2(q) + Power2(Max(0,-j + k)))*psiexp(k - nu,1))/2. + vEG2G2G2cAD*b[nu]*Power2(1/(-1 + psiexp(1,1)))*psiexp(j - k,1)*(psiexp(1 + n - q,1) + (-1 + n - q)*psiexp(1 + Max(0,-j + k),1) + Max(0,-j + k)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,-j + k),1) + (-n + q)*psiexp(2 + Max(0,-j + k),1)) + vEG2G2G2cAWAD*b[nu]*Power2(1/(-1 + psiexp(1,1)))*psiexp(j - nu,1)*(psiexp(1 + n - q,1) + (-1 + n - q)*psiexp(1 + Max(0,-j + k),1) + Max(0,-j + k)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,-j + k),1) + (-n + q)*psiexp(2 + Max(0,-j + k),1)) + vEG2G2G2cAWBD*b[nu]*Power2(1/(-1 + psiexp(1,2)))*psiexp(j - k,2)*psiexp(k - nu,1)*(psiexp(1 + n - q,2) + (-1 + n - q)*psiexp(1 + Max(0,-j + k),2) + Max(0,-j + k)*(-1 + psiexp(1,2))*psiexp(1 + Max(0,-j + k),2) + (-n + q)*psiexp(2 + Max(0,-j + k),2));
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += (azero*MMjkP1)/(n - q);
	}
}

for (int j = 1; j <= -1 + q; j++)
{
	for (int k = 1; k <= -1 + q; k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = Max(1 + j,1 + k); nu <= q; nu++)
			{
				sumInu += (j + n - nu - q)*vEG2G4c*b[nu] + (j + n - nu - q)*vEG2G4cAW*b[nu]*psiexp(-k + nu,1);
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += (azero*MMjkP1)/(n - q);
	}
}

for (int j = 1; j <= q; j++)
{
	for (int k = 1; k <= Min(-2 + j + n - q,q); k++)
	{
		MM(j, k) += (azero*vEG4G2c*b[k]*(-n + q - 2*n*q + (1 - 2*n + 2*q)*Max(0,-j + k) + Power2(n) + Power2(q) + Power2(Max(0,-j + k))))/(2.*(n - q)) + (azero*vEG4G2cAD*b[k]*Power2(1/(-1 + psiexp(1,1)))*psiexp(j - k,1)*(psiexp(1 + n - q,1) + (-1 + n - q)*psiexp(1 + Max(0,-j + k),1) + Max(0,-j + k)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,-j + k),1) + (-n + q)*psiexp(2 + Max(0,-j + k),1)))/(n - q) + (azero*vEG4G2cBD*b[k]*Power2(1/(-1 + psiexp(1,2)))*psiexp(j - k,2)*(psiexp(1 + n - q,2) + (-1 + n - q)*psiexp(1 + Max(0,-j + k),2) + Max(0,-j + k)*(-1 + psiexp(1,2))*psiexp(1 + Max(0,-j + k),2) + (-n + q)*psiexp(2 + Max(0,-j + k),2)))/(n - q);
	}
}

for (int j = 1; j <= -2 + q; j++)
{
	for (int k = Max(1,1 + j); k <= -1 + q; k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = 1 + k; nu <= q; nu++)
			{
				sumInu += (j - k + n - q)*vEG4G2c*b[nu] + (j - k + n - q)*vEG4G2cAD*b[nu]*psiexp(-k + nu,1) + (j - k + n - q)*vEG4G2cBD*b[nu]*psiexp(-k + nu,2);
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += (azero*MMjkP1)/(n - q);
	}
}

for (int j = 1; j <= -1 + q; j++)
{
	for (int k = Max(1,1 + j); k <= q; k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = 0; nu <= -1 + k; nu++)
			{
				sumInu += (j - k + n - q)*vEG2G4c*b[nu] + (j - k + n - q)*vEG2G4cAW*b[nu]*psiexp(k - nu,1);
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += (azero*MMjkP1)/(n - q);
	}
}

for (int j = 1; j <= -1 + q; j++)
{
	for (int k = Max(1,1 + j); k <= q; k++)
	{
		MM(j, k) += (azero*(j - k + n - q)*vEG6c*b[k])/(n - q);
	}
}

for (int j = 1; j <= -3 + q; j++)
{
	for (int k = Max(1,3 + j); k <= q; k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = 1 + j; nu <= Min(-2 + k,q); nu++)
			{
				sumInu += -((-1 + k - nu)*(k - 2*n - nu + 2*q)*vEG2G2G2c*b[nu])/2. + vEG2G2G2cAD*b[nu]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-k,1)*((1 - k + n + nu - q)*psiexp(1 + k,1) + (k - n - nu + q)*psiexp(2 + k,1) + (-n + q)*psiexp(2*k - nu,1) + (-1 + n - q)*psiexp(1 + 2*k - nu,1)) - ((1 - k + nu)*(-k + 2*n + nu - 2*q)*vEG2G2G2cAW*b[nu]*psiexp(-j + nu,1))/2. + vEG2G2G2cAWBD*b[nu]*Power2(1/(-1 + psiexp(1,2)))*psiexp(-k,2)*((1 - k + n + nu - q)*psiexp(1 + k,2) + (k - n - nu + q)*psiexp(2 + k,2) + (-n + q)*psiexp(2*k - nu,2) + (-1 + n - q)*psiexp(1 + 2*k - nu,2))*psiexp(-j + nu,1) + vEG2G2G2cAWAD*b[nu]*Power2(1/(-1 + psiexp(1,1)))*((1 - k + n + nu - q)*psiexp(1 + k,1) + (k - n - nu + q)*psiexp(2 + k,1) + (-n + q)*psiexp(2*k - nu,1) + (-1 + n - q)*psiexp(1 + 2*k - nu,1))*psiexp(-j - k + nu,1);
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += (azero*MMjkP1)/(n - q);
	}
}

for (int j = 1; j <= -2 + q; j++)
{
	for (int k = Max(1,2 + j); k <= q; k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = 2 + j; nu <= q; nu++)
			{
				sumInu += -(vEG2G2G2c*(j - k + 2*n - 2*q - Max(0,k - nu))*(1 + j - k + Max(0,k - nu))*b[nu])/2. - (vEG2G2G2cAWAD*(j - k + 2*n - 2*q - Max(0,k - nu))*(1 + j - k + Max(0,k - nu))*b[nu]*psiexp(-j + nu,1))/2. - vEG2G2G2cAW*b[nu]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-Max(0,k - nu),1)*((n - q)*psiexp(-j + k,1) + Max(0,k - nu)*(-1 + psiexp(1,1))*psiexp(-j + k,1) + (1 - n + q)*psiexp(1 - j + k,1) + (-1 - j + k - n + q)*psiexp(1 + Max(0,k - nu),1) + (j - k + n - q)*psiexp(2 + Max(0,k - nu),1)) + vEG2G2G2cAD*b[nu]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-k + nu,1)*((-j + k - n + q)*psiexp(-j + k,1) + (1 + j - k + n - q)*psiexp(1 - j + k,1) + (-1 + n - q)*psiexp(1 + Max(0,k - nu),1) + Max(0,k - nu)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,k - nu),1) + (-n + q)*psiexp(2 + Max(0,k - nu),1)) - vEG2G2G2cAWBD*b[nu]*Power2(1/(psiexp(1,1) - psiexp(1,2)))*psiexp(-k + nu,2)*psiexp(-Max(0,k - nu),1)*((-1 - j + k - n + q)*psiexp(1 - j + k,2)*psiexp(1 + Max(0,k - nu),1) + Max(0,k - nu)*(psiexp(1,1) - psiexp(1,2))*psiexp(-j + k,1)*psiexp(1 + Max(0,k - nu),2) + (1 - n + q)*psiexp(1 - j + k,1)*psiexp(1 + Max(0,k - nu),2) + (j - k + n - q)*psiexp(-j + k,2)*psiexp(2 + Max(0,k - nu),1) + (n - q)*psiexp(-j + k,1)*psiexp(2 + Max(0,k - nu),2));
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += (azero*MMjkP1)/(n - q);
	}
}

for (int j = 1; j <= -1 + q; j++)
{
	for (int k = 1; k <= Min(-2 + j + n - q,q); k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = 1 + j; nu <= q; nu++)
			{
				sumInu += (vEG2G2G2c*b[nu]*(-n + q - 2*n*q + (1 - 2*n + 2*q)*Max(0,-j + k) + Power2(n) + Power2(q) + Power2(Max(0,-j + k))))/2. + (vEG2G2G2cAD*b[nu]*(-n + q - 2*n*q + (1 - 2*n + 2*q)*Max(0,-j + k) + Power2(n) + Power2(q) + Power2(Max(0,-j + k)))*psiexp(-j + nu,1))/2. + vEG2G2G2cAW*b[nu]*Power2(1/(-1 + psiexp(1,1)))*psiexp(j - k,1)*(psiexp(1 + n - q,1) + (-1 + n - q)*psiexp(1 + Max(0,-j + k),1) + Max(0,-j + k)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,-j + k),1) + (-n + q)*psiexp(2 + Max(0,-j + k),1)) + vEG2G2G2cAWBD*b[nu]*Power2(1/(-1 + psiexp(1,1)))*psiexp(j - k,1)*psiexp(-j + nu,2)*(psiexp(1 + n - q,1) + (-1 + n - q)*psiexp(1 + Max(0,-j + k),1) + Max(0,-j + k)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,-j + k),1) + (-n + q)*psiexp(2 + Max(0,-j + k),1)) + vEG2G2G2cAWAD*b[nu]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-k + nu,1)*(psiexp(1 + n - q,1) + (-1 + n - q)*psiexp(1 + Max(0,-j + k),1) + Max(0,-j + k)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,-j + k),1) + (-n + q)*psiexp(2 + Max(0,-j + k),1));
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += (azero*MMjkP1)/(n - q);
	}
}

for (int j = 1; j <= -2 + q; j++)
{
	for (int k = Max(1,2 + j); k <= q; k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = 1 + j; nu <= Min(-1 + k,q); nu++)
			{
				sumInu += (-k + n + nu - q)*vEG2G4c*b[nu] + (-k + n + nu - q)*vEG2G4cAW*b[nu]*psiexp(-j + nu,1);
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += (azero*MMjkP1)/(n - q);
	}
}

for (int j = 1; j <= -1 + q; j++)
{
	for (int k = Max(1,1 + j); k <= q; k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = 1 + j; nu <= q; nu++)
			{
				sumInu += vEG4G2c*b[nu] + vEG4G2cAD*b[nu]*psiexp(-j + nu,1) + vEG4G2cBD*b[nu]*psiexp(-j + nu,2);
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += (azero*MMjkP1*(j - k + n - q))/(n - q);
	}
}

for (int j = 1; j <= -2 + q; j++)
{
	for (int k = Max(1,2 + j); k <= q; k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = 0; nu <= -1 + j; nu++)
			{
				sumInu += -((1 + j - k)*(j - k + 2*n - 2*q)*vEG2G2G2c*b[nu])/2. + vEG2G2G2cAD*b[nu]*Power2(1/(-1 + psiexp(1,1)))*(-n + q + (-1 + n - q)*psiexp(1,1) + (1 + j - k + n - q)*psiexp(1 + j - k,1) + (-j + k - n + q)*psiexp(2 + j - k,1))*psiexp(-j + k,1) - ((1 + j - k)*(j - k + 2*n - 2*q)*vEG2G2G2cAW*b[nu]*psiexp(j - nu,1))/2. + vEG2G2G2cAWBD*b[nu]*Power2(1/(-1 + psiexp(1,2)))*(-n + q + (-1 + n - q)*psiexp(1,2) + (1 + j - k + n - q)*psiexp(1 + j - k,2) + (-j + k - n + q)*psiexp(2 + j - k,2))*psiexp(-j + k,2)*psiexp(j - nu,1) + vEG2G2G2cAWAD*b[nu]*Power2(1/(-1 + psiexp(1,1)))*(-n + q + (-1 + n - q)*psiexp(1,1) + (1 + j - k + n - q)*psiexp(1 + j - k,1) + (-j + k - n + q)*psiexp(2 + j - k,1))*psiexp(k - nu,1);
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += (azero*MMjkP1)/(n - q);
	}
}

for (int j = 2; j <= q; j++)
{
	for (int k = 2; k <= q; k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = 0; nu <= Min(-2 + j,-2 + k); nu++)
			{
				sumInu += -(vEG2G2G2c*(-1 + k - nu - Max(0,-j + k))*(k - 2*n - nu + 2*q + Max(0,-j + k))*b[nu])/2. - (vEG2G2G2cAWAD*(-k + 2*n + nu - 2*q - Max(0,-j + k))*(1 - k + nu + Max(0,-j + k))*b[nu]*psiexp(j - nu,1))/2. + vEG2G2G2cAD*b[nu]*Power2(1/(-1 + psiexp(1,1)))*psiexp(j - k,1)*((k - n - nu + q)*psiexp(k - nu,1) + (1 - k + n + nu - q)*psiexp(1 + k - nu,1) + (-1 + n - q)*psiexp(1 + Max(0,-j + k),1) + Max(0,-j + k)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,-j + k),1) + (-n + q)*psiexp(2 + Max(0,-j + k),1)) + vEG2G2G2cAW*b[nu]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-Max(0,-j + k),1)*((-n + q)*psiexp(k - nu,1) - Max(0,-j + k)*(-1 + psiexp(1,1))*psiexp(k - nu,1) + (-1 + n - q)*psiexp(1 + k - nu,1) + (1 - k + n + nu - q)*psiexp(1 + Max(0,-j + k),1) + (k - n - nu + q)*psiexp(2 + Max(0,-j + k),1)) + vEG2G2G2cAWBD*b[nu]*Power2(1/(psiexp(1,1) - psiexp(1,2)))*psiexp(j - k,2)*psiexp(-Max(0,-j + k),1)*((1 - k + n + nu - q)*psiexp(1 + k - nu,2)*psiexp(1 + Max(0,-j + k),1) + Max(0,-j + k)*(-psiexp(1,1) + psiexp(1,2))*psiexp(k - nu,1)*psiexp(1 + Max(0,-j + k),2) + (-1 + n - q)*psiexp(1 + k - nu,1)*psiexp(1 + Max(0,-j + k),2) + (k - n - nu + q)*psiexp(k - nu,2)*psiexp(2 + Max(0,-j + k),1) + (-n + q)*psiexp(k - nu,1)*psiexp(2 + Max(0,-j + k),2));
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += (azero*MMjkP1)/(n - q);
	}
}

for (int j = 1; j <= q; j++)
{
	for (int k = 1; k <= Min(-3 + j + n - q,q); k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = Max(0,2 + k - n + q); nu <= -1 + j; nu++)
			{
				sumInu += (vEG2G2G2c*b[nu]*(-n + q - 2*n*q + (1 - 2*n + 2*q)*Max(0,k - nu) + Power2(n) + Power2(q) + Power2(Max(0,k - nu))))/2. + (vEG2G2G2cAD*b[nu]*(-n + q - 2*n*q + (1 - 2*n + 2*q)*Max(0,k - nu) + Power2(n) + Power2(q) + Power2(Max(0,k - nu)))*psiexp(j - nu,1))/2. + vEG2G2G2cAWAD*b[nu]*Power2(1/(-1 + psiexp(1,1)))*psiexp(j - k,1)*(psiexp(1 + n - q,1) + (-1 + n - q)*psiexp(1 + Max(0,k - nu),1) + Max(0,k - nu)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,k - nu),1) + (-n + q)*psiexp(2 + Max(0,k - nu),1)) + vEG2G2G2cAW*b[nu]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-k + nu,1)*(psiexp(1 + n - q,1) + (-1 + n - q)*psiexp(1 + Max(0,k - nu),1) + Max(0,k - nu)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,k - nu),1) + (-n + q)*psiexp(2 + Max(0,k - nu),1)) + vEG2G2G2cAWBD*b[nu]*Power2(1/(-1 + psiexp(1,1)))*psiexp(j - nu,2)*psiexp(-k + nu,1)*(psiexp(1 + n - q,1) + (-1 + n - q)*psiexp(1 + Max(0,k - nu),1) + Max(0,k - nu)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,k - nu),1) + (-n + q)*psiexp(2 + Max(0,k - nu),1));
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += (azero*MMjkP1)/(n - q);
	}
}

for (int j = 1; j <= -1 + q; j++)
{
	for (int k = Max(1,1 + j); k <= q; k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = 0; nu <= -1 + j; nu++)
			{
				sumInu += vEG2G4c*b[nu] + vEG2G4cAW*b[nu]*psiexp(j - nu,1);
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += (azero*MMjkP1*(j - k + n - q))/(n - q);
	}
}

for (int j = 1; j <= q; j++)
{
	for (int k = 1; k <= Min(-2 + j + n - q,q); k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = Max(0,1 + k - n + q); nu <= Min(-1 + j,-1 + k); nu++)
			{
				sumInu += (-k + n + nu - q)*vEG4G2c*b[nu] + (-k + n + nu - q)*vEG4G2cAD*b[nu]*psiexp(j - nu,1) + (-k + n + nu - q)*vEG4G2cBD*b[nu]*psiexp(j - nu,2);
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += (azero*MMjkP1)/(n - q);
	}
}

for (int j = 1; j <= -2 + q; j++)
{
	for (int k = Max(1,2 + j); k <= q; k++)
	{
		MM(j, k) += -(azero*(1 + j - k)*(j - k + 2*n - 2*q)*vEG4G2c*b[j])/(2.*(n - q)) + (azero*vEG4G2cAD*b[j]*Power2(1/(-1 + psiexp(1,1)))*(-n + q + (-1 + n - q)*psiexp(1,1) + (1 + j - k + n - q)*psiexp(1 + j - k,1) + (-j + k - n + q)*psiexp(2 + j - k,1))*psiexp(-j + k,1))/(n - q) + (azero*vEG4G2cBD*b[j]*Power2(1/(-1 + psiexp(1,2)))*(-n + q + (-1 + n - q)*psiexp(1,2) + (1 + j - k + n - q)*psiexp(1 + j - k,2) + (-j + k - n + q)*psiexp(2 + j - k,2))*psiexp(-j + k,2))/(n - q);
	}
}

for (int j = 1; j <= q; j++)
{
	for (int k = 1; k <= Min(-2 + j + n - q,q); k++)
	{
		MM(j, k) += (azero*vEG2G4c*b[j]*(-n + q - 2*n*q + (1 - 2*n + 2*q)*Max(0,-j + k) + Power2(n) + Power2(q) + Power2(Max(0,-j + k))))/(2.*(n - q)) + (azero*vEG2G4cAW*b[j]*Power2(1/(-1 + psiexp(1,1)))*psiexp(j - k,1)*(psiexp(1 + n - q,1) + (-1 + n - q)*psiexp(1 + Max(0,-j + k),1) + Max(0,-j + k)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,-j + k),1) + (-n + q)*psiexp(2 + Max(0,-j + k),1)))/(n - q);
	}
}

for (int j = 1; j <= -1 + q; j++)
{
	for (int k = Max(1,1 + j); k <= q; k++)
	{
		MM(j, k) += (azero*(j - k + n - q)*vEG6c*b[j])/(n - q);
	}
}

for (int j = 1; j <= -4 + q; j++)
{
	for (int k = Max(1,3 + j); k <= -1 + q; k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = 1 + j; nu <= Min(-2 + k,q); nu++)
			{
				scalar sumInuP1;
				{
					scalar sumIkappa = 0;
					for (int kappa = 1 + k; kappa <= q; kappa++)
					{
						sumIkappa += -((-1 + k - nu)*(k - 2*n - nu + 2*q)*vEG2G2G2G2c*b[kappa])/2. - ((-1 + k - nu)*(k - 2*n - nu + 2*q)*vEG2G2G2G2cAH*b[kappa]*psiexp(-k + kappa,1))/2. + vEG2G2G2G2cAD*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-k,1)*((1 - k + n + nu - q)*psiexp(1 + k,1) + (k - n - nu + q)*psiexp(2 + k,1) + (-n + q)*psiexp(2*k - nu,1) + (-1 + n - q)*psiexp(1 + 2*k - nu,1)) + vEG2G2G2G2cADAH*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-2*k + kappa,1)*((1 - k + n + nu - q)*psiexp(1 + k,1) + (k - n - nu + q)*psiexp(2 + k,1) + (-n + q)*psiexp(2*k - nu,1) + (-1 + n - q)*psiexp(1 + 2*k - nu,1)) + vEG2G2G2G2cADBH*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-k,1)*psiexp(-k + kappa,2)*((1 - k + n + nu - q)*psiexp(1 + k,1) + (k - n - nu + q)*psiexp(2 + k,1) + (-n + q)*psiexp(2*k - nu,1) + (-1 + n - q)*psiexp(1 + 2*k - nu,1)) - ((1 - k + nu)*(-k + 2*n + nu - 2*q)*vEG2G2G2G2cAW*b[kappa]*psiexp(-j + nu,1))/2. + vEG2G2G2G2cAWBD*b[kappa]*Power2(1/(-1 + psiexp(1,2)))*psiexp(-k,2)*((1 - k + n + nu - q)*psiexp(1 + k,2) + (k - n - nu + q)*psiexp(2 + k,2) + (-n + q)*psiexp(2*k - nu,2) + (-1 + n - q)*psiexp(1 + 2*k - nu,2))*psiexp(-j + nu,1) + vEG2G2G2G2cAWBDBH*b[kappa]*Power2(1/(-1 + psiexp(1,2)))*psiexp(-2*k + kappa,2)*((1 - k + n + nu - q)*psiexp(1 + k,2) + (k - n - nu + q)*psiexp(2 + k,2) + (-n + q)*psiexp(2*k - nu,2) + (-1 + n - q)*psiexp(1 + 2*k - nu,2))*psiexp(-j + nu,1) + vEG2G2G2G2cAWBDCH*b[kappa]*Power2(1/(-1 + psiexp(1,2)))*psiexp(-k,2)*psiexp(-k + kappa,3)*((1 - k + n + nu - q)*psiexp(1 + k,2) + (k - n - nu + q)*psiexp(2 + k,2) + (-n + q)*psiexp(2*k - nu,2) + (-1 + n - q)*psiexp(1 + 2*k - nu,2))*psiexp(-j + nu,1) + vEG2G2G2G2cAWAD*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*((1 - k + n + nu - q)*psiexp(1 + k,1) + (k - n - nu + q)*psiexp(2 + k,1) + (-n + q)*psiexp(2*k - nu,1) + (-1 + n - q)*psiexp(1 + 2*k - nu,1))*psiexp(-j - k + nu,1) + vEG2G2G2G2cAWADBH*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-k + kappa,2)*((1 - k + n + nu - q)*psiexp(1 + k,1) + (k - n - nu + q)*psiexp(2 + k,1) + (-n + q)*psiexp(2*k - nu,1) + (-1 + n - q)*psiexp(1 + 2*k - nu,1))*psiexp(-j - k + nu,1) + vEG2G2G2G2cAWADAH*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*((1 - k + n + nu - q)*psiexp(1 + k,1) + (k - n - nu + q)*psiexp(2 + k,1) + (-n + q)*psiexp(2*k - nu,1) + (-1 + n - q)*psiexp(1 + 2*k - nu,1))*psiexp(-j - 2*k + kappa + nu,1) - ((-1 + k - nu)*(k - 2*n - nu + 2*q)*vEG2G2G2G2cAWAH*b[kappa]*psiexp(-j - k + kappa + nu,1))/2. + vEG2G2G2G2cAWAHBD*b[kappa]*Power2(1/(-1 + psiexp(1,2)))*psiexp(-k,2)*((1 - k + n + nu - q)*psiexp(1 + k,2) + (k - n - nu + q)*psiexp(2 + k,2) + (-n + q)*psiexp(2*k - nu,2) + (-1 + n - q)*psiexp(1 + 2*k - nu,2))*psiexp(-j - k + kappa + nu,1);
					}
					sumInuP1 = sumIkappa;
				}
				sumInu += sumInuP1*b[nu];
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += MMjkP1/(n - q);
	}
}

for (int j = 1; j <= -4 + q; j++)
{
	for (int k = Max(1,4 + j); k <= q; k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = 1 + j; nu <= Min(-3 + k,q); nu++)
			{
				scalar sumInuP1;
				{
					scalar sumIkappa = 0;
					for (int kappa = Max(0,2 + nu); kappa <= -1 + k; kappa++)
					{
						sumIkappa += ((-1 + kappa - nu)*(-kappa + 2*n + nu - 2*q)*vEG2G2G2G2c*b[kappa])/2. - ((-1 + kappa - nu)*(kappa - 2*n - nu + 2*q)*vEG2G2G2G2cAH*b[kappa]*psiexp(k - kappa,1))/2. + vEG2G2G2G2cADAH*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(k - 2*kappa,1)*((1 - kappa + n + nu - q)*psiexp(1 + kappa,1) + (kappa - n - nu + q)*psiexp(2 + kappa,1) + (-n + q)*psiexp(2*kappa - nu,1) + (-1 + n - q)*psiexp(1 + 2*kappa - nu,1)) + vEG2G2G2G2cAD*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-kappa,1)*((1 - kappa + n + nu - q)*psiexp(1 + kappa,1) + (kappa - n - nu + q)*psiexp(2 + kappa,1) + (-n + q)*psiexp(2*kappa - nu,1) + (-1 + n - q)*psiexp(1 + 2*kappa - nu,1)) + vEG2G2G2G2cADBH*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(k - kappa,2)*psiexp(-kappa,1)*((1 - kappa + n + nu - q)*psiexp(1 + kappa,1) + (kappa - n - nu + q)*psiexp(2 + kappa,1) + (-n + q)*psiexp(2*kappa - nu,1) + (-1 + n - q)*psiexp(1 + 2*kappa - nu,1)) - ((1 - kappa + nu)*(-kappa + 2*n + nu - 2*q)*vEG2G2G2G2cAW*b[kappa]*psiexp(-j + nu,1))/2. + vEG2G2G2G2cAWBDBH*b[kappa]*Power2(1/(-1 + psiexp(1,2)))*psiexp(k - 2*kappa,2)*((1 - kappa + n + nu - q)*psiexp(1 + kappa,2) + (kappa - n - nu + q)*psiexp(2 + kappa,2) + (-n + q)*psiexp(2*kappa - nu,2) + (-1 + n - q)*psiexp(1 + 2*kappa - nu,2))*psiexp(-j + nu,1) + vEG2G2G2G2cAWBD*b[kappa]*Power2(1/(-1 + psiexp(1,2)))*psiexp(-kappa,2)*((1 - kappa + n + nu - q)*psiexp(1 + kappa,2) + (kappa - n - nu + q)*psiexp(2 + kappa,2) + (-n + q)*psiexp(2*kappa - nu,2) + (-1 + n - q)*psiexp(1 + 2*kappa - nu,2))*psiexp(-j + nu,1) + vEG2G2G2G2cAWBDCH*b[kappa]*Power2(1/(-1 + psiexp(1,2)))*psiexp(k - kappa,3)*psiexp(-kappa,2)*((1 - kappa + n + nu - q)*psiexp(1 + kappa,2) + (kappa - n - nu + q)*psiexp(2 + kappa,2) + (-n + q)*psiexp(2*kappa - nu,2) + (-1 + n - q)*psiexp(1 + 2*kappa - nu,2))*psiexp(-j + nu,1) + vEG2G2G2G2cAWADAH*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*((1 - kappa + n + nu - q)*psiexp(1 + kappa,1) + (kappa - n - nu + q)*psiexp(2 + kappa,1) + (-n + q)*psiexp(2*kappa - nu,1) + (-1 + n - q)*psiexp(1 + 2*kappa - nu,1))*psiexp(-j + k - 2*kappa + nu,1) + vEG2G2G2G2cAWAD*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*((1 - kappa + n + nu - q)*psiexp(1 + kappa,1) + (kappa - n - nu + q)*psiexp(2 + kappa,1) + (-n + q)*psiexp(2*kappa - nu,1) + (-1 + n - q)*psiexp(1 + 2*kappa - nu,1))*psiexp(-j - kappa + nu,1) + vEG2G2G2G2cAWADBH*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(k - kappa,2)*((1 - kappa + n + nu - q)*psiexp(1 + kappa,1) + (kappa - n - nu + q)*psiexp(2 + kappa,1) + (-n + q)*psiexp(2*kappa - nu,1) + (-1 + n - q)*psiexp(1 + 2*kappa - nu,1))*psiexp(-j - kappa + nu,1) - ((-1 + kappa - nu)*(kappa - 2*n - nu + 2*q)*vEG2G2G2G2cAWAH*b[kappa]*psiexp(-j + k - kappa + nu,1))/2. + vEG2G2G2G2cAWAHBD*b[kappa]*Power2(1/(-1 + psiexp(1,2)))*psiexp(-kappa,2)*((1 - kappa + n + nu - q)*psiexp(1 + kappa,2) + (kappa - n - nu + q)*psiexp(2 + kappa,2) + (-n + q)*psiexp(2*kappa - nu,2) + (-1 + n - q)*psiexp(1 + 2*kappa - nu,2))*psiexp(-j + k - kappa + nu,1);
					}
					sumInuP1 = sumIkappa;
				}
				sumInu += sumInuP1*b[nu];
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += MMjkP1/(n - q);
	}
}

for (int j = 1; j <= -4 + q; j++)
{
	for (int k = Max(2,4 + j); k <= q; k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = 2 + j; nu <= Min3(-2 + k,j + k,q); nu++)
			{
				scalar sumInuP1;
				{
					scalar sumIkappa = 0;
					for (int kappa = Max(0,2 + j); kappa <= j + k - nu; kappa++)
					{
						sumIkappa += -(vEG2G2G2G2c*(j - kappa + 2*n - 2*q - Max(0,kappa - nu))*(1 + j - kappa + Max(0,kappa - nu))*b[kappa])/2. - (vEG2G2G2G2cADAH*(-1 - j + kappa - Max(0,kappa - nu))*(-j + kappa - 2*n + 2*q + Max(0,kappa - nu))*b[kappa]*psiexp(k - kappa,1))/2. - (vEG2G2G2G2cAWAD*(j - kappa + 2*n - 2*q - Max(0,kappa - nu))*(1 + j - kappa + Max(0,kappa - nu))*b[kappa]*psiexp(-j + nu,1))/2. - vEG2G2G2G2cAWADAH*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(k - kappa - Max(0,kappa - nu),1)*((n - q)*psiexp(-j + kappa,1) + Max(0,kappa - nu)*(-1 + psiexp(1,1))*psiexp(-j + kappa,1) + (1 - n + q)*psiexp(1 - j + kappa,1) + (-1 - j + kappa - n + q)*psiexp(1 + Max(0,kappa - nu),1) + (j - kappa + n - q)*psiexp(2 + Max(0,kappa - nu),1)) - vEG2G2G2G2cAH*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(j + k - kappa - nu - Max(0,kappa - nu),1)*((n - q)*psiexp(-j + kappa,1) + Max(0,kappa - nu)*(-1 + psiexp(1,1))*psiexp(-j + kappa,1) + (1 - n + q)*psiexp(1 - j + kappa,1) + (-1 - j + kappa - n + q)*psiexp(1 + Max(0,kappa - nu),1) + (j - kappa + n - q)*psiexp(2 + Max(0,kappa - nu),1)) - vEG2G2G2G2cAW*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-Max(0,kappa - nu),1)*((n - q)*psiexp(-j + kappa,1) + Max(0,kappa - nu)*(-1 + psiexp(1,1))*psiexp(-j + kappa,1) + (1 - n + q)*psiexp(1 - j + kappa,1) + (-1 - j + kappa - n + q)*psiexp(1 + Max(0,kappa - nu),1) + (j - kappa + n - q)*psiexp(2 + Max(0,kappa - nu),1)) + vEG2G2G2G2cAD*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-kappa + nu,1)*((-j + kappa - n + q)*psiexp(-j + kappa,1) + (1 + j - kappa + n - q)*psiexp(1 - j + kappa,1) + (-1 + n - q)*psiexp(1 + Max(0,kappa - nu),1) + Max(0,kappa - nu)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,kappa - nu),1) + (-n + q)*psiexp(2 + Max(0,kappa - nu),1)) + vEG2G2G2G2cAWBDBH*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(k - kappa,2)*psiexp(-Max(0,kappa - nu),1)*((-n + q)*psiexp(-j + kappa,1) - Max(0,kappa - nu)*(-1 + psiexp(1,1))*psiexp(-j + kappa,1) + (-1 + n - q)*psiexp(1 - j + kappa,1) + (1 + j - kappa + n - q)*psiexp(1 + Max(0,kappa - nu),1) + (-j + kappa - n + q)*psiexp(2 + Max(0,kappa - nu),1)) - vEG2G2G2G2cAWADBH*b[kappa]*Power2(1/(-1 + psiexp(1,2)))*psiexp(-j + nu,1)*psiexp(j + k - kappa - nu - Max(0,kappa - nu),2)*((n - q)*psiexp(-j + kappa,2) + Max(0,kappa - nu)*(-1 + psiexp(1,2))*psiexp(-j + kappa,2) + (1 - n + q)*psiexp(1 - j + kappa,2) + (-1 - j + kappa - n + q)*psiexp(1 + Max(0,kappa - nu),2) + (j - kappa + n - q)*psiexp(2 + Max(0,kappa - nu),2)) - vEG2G2G2G2cAWBD*b[kappa]*Power2(1/(psiexp(1,1) - psiexp(1,2)))*psiexp(-kappa + nu,2)*psiexp(-Max(0,kappa - nu),1)*((-1 - j + kappa - n + q)*psiexp(1 - j + kappa,2)*psiexp(1 + Max(0,kappa - nu),1) + Max(0,kappa - nu)*(psiexp(1,1) - psiexp(1,2))*psiexp(-j + kappa,1)*psiexp(1 + Max(0,kappa - nu),2) + (1 - n + q)*psiexp(1 - j + kappa,1)*psiexp(1 + Max(0,kappa - nu),2) + (j - kappa + n - q)*psiexp(-j + kappa,2)*psiexp(2 + Max(0,kappa - nu),1) + (n - q)*psiexp(-j + kappa,1)*psiexp(2 + Max(0,kappa - nu),2)) + vEG2G2G2G2cADBH*b[kappa]*Power2(1/(psiexp(1,1) - psiexp(1,2)))*psiexp(-kappa + nu,1)*psiexp(j + k - kappa - nu - Max(0,kappa - nu),2)*(Max(0,kappa - nu)*(psiexp(1,1) - psiexp(1,2))*psiexp(-j + kappa,2)*psiexp(1 + Max(0,kappa - nu),1) + (-1 + n - q)*psiexp(1 - j + kappa,2)*psiexp(1 + Max(0,kappa - nu),1) + (1 + j - kappa + n - q)*psiexp(1 - j + kappa,1)*psiexp(1 + Max(0,kappa - nu),2) + (-n + q)*psiexp(-j + kappa,2)*psiexp(2 + Max(0,kappa - nu),1) + (-j + kappa - n + q)*psiexp(-j + kappa,1)*psiexp(2 + Max(0,kappa - nu),2)) + vEG2G2G2G2cAWAH*b[kappa]*Power2(1/(-1 + psiexp(2,1)))*psiexp(j + k - kappa - nu - 2*Max(0,kappa - nu),1)*((-n + q)*psiexp(2*(-j + kappa),1) - Max(0,kappa - nu)*(-1 + psiexp(2,1))*psiexp(2*(-j + kappa),1) + (-1 + n - q)*psiexp(2*(1 - j + kappa),1) + (1 + j - kappa + n - q)*psiexp(2*(1 + Max(0,kappa - nu)),1) + (-j + kappa - n + q)*psiexp(2*(2 + Max(0,kappa - nu)),1)) - vEG2G2G2G2cAWBDCH*b[kappa]*Power2(1/(psiexp(1,2) - psiexp(1,1)*psiexp(1,3)))*psiexp(-j - kappa + nu,2)*psiexp(-j - Max(0,kappa - nu),1)*psiexp(k - kappa - nu - Max(0,kappa - nu),3)*(Max(0,kappa - nu)*(-psiexp(1,2) + psiexp(1,1)*psiexp(1,3))*psiexp(kappa,1)*psiexp(kappa,3)*psiexp(1 + j + Max(0,kappa - nu),2) + (1 - n + q)*psiexp(1 + kappa,1)*psiexp(1 + kappa,3)*psiexp(1 + j + Max(0,kappa - nu),2) + (-1 - j + kappa - n + q)*psiexp(1 + kappa,2)*psiexp(1 + j + Max(0,kappa - nu),1)*psiexp(1 + j + Max(0,kappa - nu),3) + (n - q)*psiexp(kappa,1)*psiexp(kappa,3)*psiexp(2 + j + Max(0,kappa - nu),2) + (j - kappa + n - q)*psiexp(kappa,2)*psiexp(2 + j + Max(0,kappa - nu),1)*psiexp(2 + j + Max(0,kappa - nu),3)) - vEG2G2G2G2cAWAHBD*b[kappa]*Power2(1/(-psiexp(1,2) + psiexp(2,1)))*psiexp(-kappa + nu,2)*psiexp(j + k - kappa - nu - 2*Max(0,kappa - nu),1)*(Max(0,kappa - nu)*(-psiexp(1,2) + psiexp(2,1))*psiexp(2*(-j + kappa),1)*psiexp(1 + Max(0,kappa - nu),2) + (1 - n + q)*psiexp(2 - 2*j + 2*kappa,1)*psiexp(1 + Max(0,kappa - nu),2) + (-1 - j + kappa - n + q)*psiexp(1 - j + kappa,2)*psiexp(2*(1 + Max(0,kappa - nu)),1) + (n - q)*psiexp(2*(-j + kappa),1)*psiexp(2 + Max(0,kappa - nu),2) + (j - kappa + n - q)*psiexp(-j + kappa,2)*psiexp(4 + 2*Max(0,kappa - nu),1));
					}
					sumInuP1 = sumIkappa;
				}
				sumInu += sumInuP1*b[nu];
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += MMjkP1/(n - q);
	}
}

for (int j = 1; j <= -5 + q; j++)
{
	for (int k = Max(2,5 + j); k <= q; k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = 3 + j; nu <= Min(-2 + k,q); nu++)
			{
				scalar sumInuP1;
				{
					scalar sumIkappa = 0;
					for (int kappa = Max(0,1 + j + k - nu); kappa <= -2 + k; kappa++)
					{
						sumIkappa += -(vEG2G2G2G2c*(-1 + k - nu - Max(0,kappa - nu))*(k - 2*n - nu + 2*q + Max(0,kappa - nu))*b[kappa])/2. - (vEG2G2G2G2cADAH*(-1 + k - nu - Max(0,kappa - nu))*(k - 2*n - nu + 2*q + Max(0,kappa - nu))*b[kappa]*psiexp(k - kappa,1))/2. - (vEG2G2G2G2cAWAD*(-k + 2*n + nu - 2*q - Max(0,kappa - nu))*(1 - k + nu + Max(0,kappa - nu))*b[kappa]*psiexp(-j + nu,1))/2. + vEG2G2G2G2cAD*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-kappa + nu,1)*((k - n - nu + q)*psiexp(k - nu,1) + (1 - k + n + nu - q)*psiexp(1 + k - nu,1) + (-1 + n - q)*psiexp(1 + Max(0,kappa - nu),1) + Max(0,kappa - nu)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,kappa - nu),1) + (-n + q)*psiexp(2 + Max(0,kappa - nu),1)) + vEG2G2G2G2cAWADAH*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-j + nu - Max(0,kappa - nu),1)*((-n + q)*psiexp(k - nu,1) - Max(0,kappa - nu)*(-1 + psiexp(1,1))*psiexp(k - nu,1) + (-1 + n - q)*psiexp(1 + k - nu,1) + (1 - k + n + nu - q)*psiexp(1 + Max(0,kappa - nu),1) + (k - n - nu + q)*psiexp(2 + Max(0,kappa - nu),1)) + vEG2G2G2G2cAW*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-j - k + kappa + nu - Max(0,kappa - nu),1)*((-n + q)*psiexp(k - nu,1) - Max(0,kappa - nu)*(-1 + psiexp(1,1))*psiexp(k - nu,1) + (-1 + n - q)*psiexp(1 + k - nu,1) + (1 - k + n + nu - q)*psiexp(1 + Max(0,kappa - nu),1) + (k - n - nu + q)*psiexp(2 + Max(0,kappa - nu),1)) + vEG2G2G2G2cAWBDBH*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(k - kappa,2)*psiexp(-j - k + kappa + nu - Max(0,kappa - nu),1)*((-n + q)*psiexp(k - nu,1) - Max(0,kappa - nu)*(-1 + psiexp(1,1))*psiexp(k - nu,1) + (-1 + n - q)*psiexp(1 + k - nu,1) + (1 - k + n + nu - q)*psiexp(1 + Max(0,kappa - nu),1) + (k - n - nu + q)*psiexp(2 + Max(0,kappa - nu),1)) + vEG2G2G2G2cAH*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-Max(0,kappa - nu),1)*((-n + q)*psiexp(k - nu,1) - Max(0,kappa - nu)*(-1 + psiexp(1,1))*psiexp(k - nu,1) + (-1 + n - q)*psiexp(1 + k - nu,1) + (1 - k + n + nu - q)*psiexp(1 + Max(0,kappa - nu),1) + (k - n - nu + q)*psiexp(2 + Max(0,kappa - nu),1)) - vEG2G2G2G2cAWADBH*b[kappa]*Power2(1/(-1 + psiexp(1,2)))*psiexp(-j + nu,1)*psiexp(-Max(0,kappa - nu),2)*((n - q)*psiexp(k - nu,2) + Max(0,kappa - nu)*(-1 + psiexp(1,2))*psiexp(k - nu,2) + (1 - n + q)*psiexp(1 + k - nu,2) + (-1 + k - n - nu + q)*psiexp(1 + Max(0,kappa - nu),2) + (-k + n + nu - q)*psiexp(2 + Max(0,kappa - nu),2)) - vEG2G2G2G2cADBH*b[kappa]*Power2(1/(psiexp(1,1) - psiexp(1,2)))*psiexp(-kappa + nu,1)*psiexp(-Max(0,kappa - nu),2)*(-(Max(0,kappa - nu)*(psiexp(1,1) - psiexp(1,2))*psiexp(k - nu,2)*psiexp(1 + Max(0,kappa - nu),1)) + (1 - n + q)*psiexp(1 + k - nu,2)*psiexp(1 + Max(0,kappa - nu),1) + (-1 + k - n - nu + q)*psiexp(1 + k - nu,1)*psiexp(1 + Max(0,kappa - nu),2) + (n - q)*psiexp(k - nu,2)*psiexp(2 + Max(0,kappa - nu),1) + (-k + n + nu - q)*psiexp(k - nu,1)*psiexp(2 + Max(0,kappa - nu),2)) + vEG2G2G2G2cAWBD*b[kappa]*Power2(1/(psiexp(1,1) - psiexp(1,2)))*psiexp(-kappa + nu,2)*psiexp(-j - k + kappa + nu - Max(0,kappa - nu),1)*((1 - k + n + nu - q)*psiexp(1 + k - nu,2)*psiexp(1 + Max(0,kappa - nu),1) + Max(0,kappa - nu)*(-psiexp(1,1) + psiexp(1,2))*psiexp(k - nu,1)*psiexp(1 + Max(0,kappa - nu),2) + (-1 + n - q)*psiexp(1 + k - nu,1)*psiexp(1 + Max(0,kappa - nu),2) + (k - n - nu + q)*psiexp(k - nu,2)*psiexp(2 + Max(0,kappa - nu),1) + (-n + q)*psiexp(k - nu,1)*psiexp(2 + Max(0,kappa - nu),2)) + vEG2G2G2G2cAWAH*b[kappa]*Power2(1/(-1 + psiexp(2,1)))*psiexp(-j - k + kappa + nu - 2*Max(0,kappa - nu),1)*((-n + q)*psiexp(2*(k - nu),1) - Max(0,kappa - nu)*(-1 + psiexp(2,1))*psiexp(2*(k - nu),1) + (-1 + n - q)*psiexp(2*(1 + k - nu),1) + (1 - k + n + nu - q)*psiexp(2*(1 + Max(0,kappa - nu)),1) + (k - n - nu + q)*psiexp(2*(2 + Max(0,kappa - nu)),1)) + vEG2G2G2G2cAWBDCH*b[kappa]*Power2(1/(psiexp(1,2) - psiexp(1,1)*psiexp(1,3)))*psiexp(-kappa,2)*psiexp(-j - k + kappa - Max(0,kappa - nu),1)*psiexp(-nu - Max(0,kappa - nu),3)*(Max(0,kappa - nu)*(psiexp(1,2) - psiexp(1,1)*psiexp(1,3))*psiexp(k,1)*psiexp(k,3)*psiexp(1 + nu + Max(0,kappa - nu),2) + (-1 + n - q)*psiexp(1 + k,1)*psiexp(1 + k,3)*psiexp(1 + nu + Max(0,kappa - nu),2) + (1 - k + n + nu - q)*psiexp(1 + k,2)*psiexp(1 + nu + Max(0,kappa - nu),1)*psiexp(1 + nu + Max(0,kappa - nu),3) + (-n + q)*psiexp(k,1)*psiexp(k,3)*psiexp(2 + nu + Max(0,kappa - nu),2) + (k - n - nu + q)*psiexp(k,2)*psiexp(2 + nu + Max(0,kappa - nu),1)*psiexp(2 + nu + Max(0,kappa - nu),3)) + vEG2G2G2G2cAWAHBD*b[kappa]*Power2(1/(-psiexp(1,2) + psiexp(2,1)))*psiexp(-kappa + nu,2)*psiexp(-j - k + kappa + nu - 2*Max(0,kappa - nu),1)*(Max(0,kappa - nu)*(psiexp(1,2) - psiexp(2,1))*psiexp(2*k - 2*nu,1)*psiexp(1 + Max(0,kappa - nu),2) + (-1 + n - q)*psiexp(2*(1 + k - nu),1)*psiexp(1 + Max(0,kappa - nu),2) + (1 - k + n + nu - q)*psiexp(1 + k - nu,2)*psiexp(2*(1 + Max(0,kappa - nu)),1) + (-n + q)*psiexp(2*(k - nu),1)*psiexp(2 + Max(0,kappa - nu),2) + (k - n - nu + q)*psiexp(k - nu,2)*psiexp(4 + 2*Max(0,kappa - nu),1));
					}
					sumInuP1 = sumIkappa;
				}
				sumInu += sumInuP1*b[nu];
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += MMjkP1/(n - q);
	}
}

for (int j = 1; j <= -3 + q; j++)
{
	for (int k = Max(3,3 + j); k <= q; k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = 1 + j; nu <= Min3(-2 + k,-2 + j + k,q); nu++)
			{
				scalar sumInuP1;
				{
					scalar sumIkappa = 0;
					for (int kappa = 0; kappa <= -2 + j + k - nu; kappa++)
					{
						sumIkappa += -(vEG2G2G2G2c*(-1 + k - nu - Max(0,-j + kappa))*(k - 2*n - nu + 2*q + Max(0,-j + kappa))*b[kappa])/2. - (vEG2G2G2G2cAWADAH*(-1 + k - nu - Max(0,-j + kappa))*(k - 2*n - nu + 2*q + Max(0,-j + kappa))*b[kappa]*psiexp(k - kappa,1))/2. - (vEG2G2G2G2cAWAH*(-1 + k - nu - Max(0,-j + kappa))*(k - 2*n - nu + 2*q + Max(0,-j + kappa))*b[kappa]*psiexp(j + k - kappa - nu,1))/2. - (vEG2G2G2G2cAD*(-k + 2*n + nu - 2*q - Max(0,-j + kappa))*(1 - k + nu + Max(0,-j + kappa))*b[kappa]*psiexp(-j + nu,1))/2. - (vEG2G2G2G2cAWAHBD*(-1 + k - nu - Max(0,-j + kappa))*(k - 2*n - nu + 2*q + Max(0,-j + kappa))*b[kappa]*psiexp(j + k - kappa - nu,1)*psiexp(-j + nu,2))/2. + vEG2G2G2G2cAW*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(j - kappa,1)*((k - n - nu + q)*psiexp(k - nu,1) + (1 - k + n + nu - q)*psiexp(1 + k - nu,1) + (-1 + n - q)*psiexp(1 + Max(0,-j + kappa),1) + Max(0,-j + kappa)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,-j + kappa),1) + (-n + q)*psiexp(2 + Max(0,-j + kappa),1)) + vEG2G2G2G2cAWBD*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(j - kappa,1)*psiexp(-j + nu,2)*((k - n - nu + q)*psiexp(k - nu,1) + (1 - k + n + nu - q)*psiexp(1 + k - nu,1) + (-1 + n - q)*psiexp(1 + Max(0,-j + kappa),1) + Max(0,-j + kappa)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,-j + kappa),1) + (-n + q)*psiexp(2 + Max(0,-j + kappa),1)) + vEG2G2G2G2cAWAD*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-kappa + nu,1)*((k - n - nu + q)*psiexp(k - nu,1) + (1 - k + n + nu - q)*psiexp(1 + k - nu,1) + (-1 + n - q)*psiexp(1 + Max(0,-j + kappa),1) + Max(0,-j + kappa)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,-j + kappa),1) + (-n + q)*psiexp(2 + Max(0,-j + kappa),1)) + vEG2G2G2G2cADAH*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-j + nu - Max(0,-j + kappa),1)*((-n + q)*psiexp(k - nu,1) - Max(0,-j + kappa)*(-1 + psiexp(1,1))*psiexp(k - nu,1) + (-1 + n - q)*psiexp(1 + k - nu,1) + (1 - k + n + nu - q)*psiexp(1 + Max(0,-j + kappa),1) + (k - n - nu + q)*psiexp(2 + Max(0,-j + kappa),1)) + vEG2G2G2G2cAH*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-Max(0,-j + kappa),1)*((-n + q)*psiexp(k - nu,1) - Max(0,-j + kappa)*(-1 + psiexp(1,1))*psiexp(k - nu,1) + (-1 + n - q)*psiexp(1 + k - nu,1) + (1 - k + n + nu - q)*psiexp(1 + Max(0,-j + kappa),1) + (k - n - nu + q)*psiexp(2 + Max(0,-j + kappa),1)) - vEG2G2G2G2cADBH*b[kappa]*Power2(1/(-1 + psiexp(1,2)))*psiexp(-j + nu,1)*psiexp(-Max(0,-j + kappa),2)*((n - q)*psiexp(k - nu,2) + Max(0,-j + kappa)*(-1 + psiexp(1,2))*psiexp(k - nu,2) + (1 - n + q)*psiexp(1 + k - nu,2) + (-1 + k - n - nu + q)*psiexp(1 + Max(0,-j + kappa),2) + (-k + n + nu - q)*psiexp(2 + Max(0,-j + kappa),2)) - vEG2G2G2G2cAWADBH*b[kappa]*Power2(1/(psiexp(1,1) - psiexp(1,2)))*psiexp(-kappa + nu,1)*psiexp(-Max(0,-j + kappa),2)*(-(Max(0,-j + kappa)*(psiexp(1,1) - psiexp(1,2))*psiexp(k - nu,2)*psiexp(1 + Max(0,-j + kappa),1)) + (1 - n + q)*psiexp(1 + k - nu,2)*psiexp(1 + Max(0,-j + kappa),1) + (-1 + k - n - nu + q)*psiexp(1 + k - nu,1)*psiexp(1 + Max(0,-j + kappa),2) + (n - q)*psiexp(k - nu,2)*psiexp(2 + Max(0,-j + kappa),1) + (-k + n + nu - q)*psiexp(k - nu,1)*psiexp(2 + Max(0,-j + kappa),2)) + vEG2G2G2G2cAWBDBH*b[kappa]*Power2(1/(psiexp(1,1) - psiexp(1,2)))*psiexp(j - kappa,1)*psiexp(-j + nu - Max(0,-j + kappa),2)*(Max(0,-j + kappa)*(psiexp(1,1) - psiexp(1,2))*psiexp(k - nu,2)*psiexp(1 + Max(0,-j + kappa),1) + (-1 + n - q)*psiexp(1 + k - nu,2)*psiexp(1 + Max(0,-j + kappa),1) + (1 - k + n + nu - q)*psiexp(1 + k - nu,1)*psiexp(1 + Max(0,-j + kappa),2) + (-n + q)*psiexp(k - nu,2)*psiexp(2 + Max(0,-j + kappa),1) + (k - n - nu + q)*psiexp(k - nu,1)*psiexp(2 + Max(0,-j + kappa),2)) - vEG2G2G2G2cAWBDCH*b[kappa]*Power2(1/(psiexp(1,1) - psiexp(1,3)))*psiexp(j - kappa,1)*psiexp(-j + nu,2)*psiexp(-Max(0,-j + kappa),3)*(-(Max(0,-j + kappa)*(psiexp(1,1) - psiexp(1,3))*psiexp(k - nu,3)*psiexp(1 + Max(0,-j + kappa),1)) + (1 - n + q)*psiexp(1 + k - nu,3)*psiexp(1 + Max(0,-j + kappa),1) + (-1 + k - n - nu + q)*psiexp(1 + k - nu,1)*psiexp(1 + Max(0,-j + kappa),3) + (n - q)*psiexp(k - nu,3)*psiexp(2 + Max(0,-j + kappa),1) + (-k + n + nu - q)*psiexp(k - nu,1)*psiexp(2 + Max(0,-j + kappa),3));
					}
					sumInuP1 = sumIkappa;
				}
				sumInu += sumInuP1*b[nu];
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += MMjkP1/(n - q);
	}
}

for (int j = 1; j <= -3 + q; j++)
{
	for (int k = Max(1,3 + j); k <= q; k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = 1 + j; nu <= Min(-2 + k,q); nu++)
			{
				sumInu += -((-1 + k - nu)*(k - 2*n - nu + 2*q)*vEG2G2G4c*b[nu])/2. + vEG2G2G4cAD*b[nu]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-k,1)*((1 - k + n + nu - q)*psiexp(1 + k,1) + (k - n - nu + q)*psiexp(2 + k,1) + (-n + q)*psiexp(2*k - nu,1) + (-1 + n - q)*psiexp(1 + 2*k - nu,1)) - ((1 - k + nu)*(-k + 2*n + nu - 2*q)*vEG2G2G4cAW*b[nu]*psiexp(-j + nu,1))/2. + vEG2G2G4cAWBD*b[nu]*Power2(1/(-1 + psiexp(1,2)))*psiexp(-k,2)*((1 - k + n + nu - q)*psiexp(1 + k,2) + (k - n - nu + q)*psiexp(2 + k,2) + (-n + q)*psiexp(2*k - nu,2) + (-1 + n - q)*psiexp(1 + 2*k - nu,2))*psiexp(-j + nu,1) + vEG2G2G4cAWAD*b[nu]*Power2(1/(-1 + psiexp(1,1)))*((1 - k + n + nu - q)*psiexp(1 + k,1) + (k - n - nu + q)*psiexp(2 + k,1) + (-n + q)*psiexp(2*k - nu,1) + (-1 + n - q)*psiexp(1 + 2*k - nu,1))*psiexp(-j - k + nu,1);
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += (MMjkP1*b[k])/(n - q);
	}
}

for (int j = 1; j <= -3 + q; j++)
{
	for (int k = Max(1,3 + j); k <= q; k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = 1 + j; nu <= Min(-2 + k,q); nu++)
			{
				scalar sumInuP1;
				{
					scalar sumIkappa = 0;
					for (int kappa = Max(0,1 + nu); kappa <= -1 + k; kappa++)
					{
						sumIkappa += (-kappa + n + nu - q)*vEG2G4G2c*b[kappa] - (kappa - n - nu + q)*vEG2G4G2cAH*b[kappa]*psiexp(k - kappa,1) - (kappa - n - nu + q)*vEG2G4G2cBH*b[kappa]*psiexp(k - kappa,2) + (-kappa + n + nu - q)*vEG2G4G2cAW*b[kappa]*psiexp(-j + nu,1) + (-kappa + n + nu - q)*vEG2G4G2cAWBH*b[kappa]*psiexp(k - kappa,2)*psiexp(-j + nu,1) + (-kappa + n + nu - q)*vEG2G4G2cAWCH*b[kappa]*psiexp(k - kappa,3)*psiexp(-j + nu,1) - (kappa - n - nu + q)*vEG2G4G2cAWAH*b[kappa]*psiexp(-j + k - kappa + nu,1);
					}
					sumInuP1 = sumIkappa;
				}
				sumInu += sumInuP1*b[nu];
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += MMjkP1/(n - q);
	}
}

for (int j = 1; j <= -3 + q; j++)
{
	for (int k = Max(2,3 + j); k <= q; k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = 1 + j; nu <= Min3(-2 + k,-1 + j + k,q); nu++)
			{
				scalar sumInuP1;
				{
					scalar sumIkappa = 0;
					for (int kappa = Max(0,1 + j); kappa <= -1 + j + k - nu; kappa++)
					{
						sumIkappa += (j - kappa + n - q)*vEG4G2G2c*b[kappa] + (j - kappa + n - q)*vEG4G2G2cADAH*b[kappa]*psiexp(k - kappa,1) + (j - kappa + n - q)*vEG4G2G2cBDBH*b[kappa]*psiexp(k - kappa,2) + (j - kappa + n - q)*vEG4G2G2cAH*b[kappa]*psiexp(j + k - kappa - nu,1) + (j - kappa + n - q)*vEG4G2G2cAD*b[kappa]*psiexp(-j + nu,1) + (j - kappa + n - q)*vEG4G2G2cADBH*b[kappa]*psiexp(j + k - kappa - nu,2)*psiexp(-j + nu,1) + (j - kappa + n - q)*vEG4G2G2cBD*b[kappa]*psiexp(-j + nu,2) + (j - kappa + n - q)*vEG4G2G2cAHBD*b[kappa]*psiexp(j + k - kappa - nu,1)*psiexp(-j + nu,2) + (j - kappa + n - q)*vEG4G2G2cBDCH*b[kappa]*psiexp(j + k - kappa - nu,3)*psiexp(-j + nu,2);
					}
					sumInuP1 = sumIkappa;
				}
				sumInu += sumInuP1*b[nu];
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += MMjkP1/(n - q);
	}
}

for (int j = 1; j <= -5 + q; j++)
{
	for (int k = Max(1,2 + j); k <= -3 + q; k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = 2 + j; nu <= Min(q,-1 + j - k + q); nu++)
			{
				scalar sumInuP1;
				{
					scalar sumIkappa = 0;
					for (int kappa = 1 - j + k + nu; kappa <= q; kappa++)
					{
						sumIkappa += -(vEG2G2G2G2c*(j - k + 2*n - 2*q - Max(0,k - nu))*(1 + j - k + Max(0,k - nu))*b[kappa])/2. - (vEG2G2G2G2cADAH*(j - k + 2*n - 2*q - Max(0,k - nu))*(1 + j - k + Max(0,k - nu))*b[kappa]*psiexp(-k + kappa,1))/2. - (vEG2G2G2G2cAWAD*(j - k + 2*n - 2*q - Max(0,k - nu))*(1 + j - k + Max(0,k - nu))*b[kappa]*psiexp(-j + nu,1))/2. - vEG2G2G2G2cAWADAH*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-k + kappa - Max(0,k - nu),1)*((n - q)*psiexp(-j + k,1) + Max(0,k - nu)*(-1 + psiexp(1,1))*psiexp(-j + k,1) + (1 - n + q)*psiexp(1 - j + k,1) + (-1 - j + k - n + q)*psiexp(1 + Max(0,k - nu),1) + (j - k + n - q)*psiexp(2 + Max(0,k - nu),1)) - vEG2G2G2G2cAH*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(j - k + kappa - nu - Max(0,k - nu),1)*((n - q)*psiexp(-j + k,1) + Max(0,k - nu)*(-1 + psiexp(1,1))*psiexp(-j + k,1) + (1 - n + q)*psiexp(1 - j + k,1) + (-1 - j + k - n + q)*psiexp(1 + Max(0,k - nu),1) + (j - k + n - q)*psiexp(2 + Max(0,k - nu),1)) - vEG2G2G2G2cAW*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-Max(0,k - nu),1)*((n - q)*psiexp(-j + k,1) + Max(0,k - nu)*(-1 + psiexp(1,1))*psiexp(-j + k,1) + (1 - n + q)*psiexp(1 - j + k,1) + (-1 - j + k - n + q)*psiexp(1 + Max(0,k - nu),1) + (j - k + n - q)*psiexp(2 + Max(0,k - nu),1)) + vEG2G2G2G2cAD*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-k + nu,1)*((-j + k - n + q)*psiexp(-j + k,1) + (1 + j - k + n - q)*psiexp(1 - j + k,1) + (-1 + n - q)*psiexp(1 + Max(0,k - nu),1) + Max(0,k - nu)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,k - nu),1) + (-n + q)*psiexp(2 + Max(0,k - nu),1)) + vEG2G2G2G2cAWBDBH*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-k + kappa,2)*psiexp(-Max(0,k - nu),1)*((-n + q)*psiexp(-j + k,1) - Max(0,k - nu)*(-1 + psiexp(1,1))*psiexp(-j + k,1) + (-1 + n - q)*psiexp(1 - j + k,1) + (1 + j - k + n - q)*psiexp(1 + Max(0,k - nu),1) + (-j + k - n + q)*psiexp(2 + Max(0,k - nu),1)) - vEG2G2G2G2cAWADBH*b[kappa]*Power2(1/(-1 + psiexp(1,2)))*psiexp(-j + nu,1)*psiexp(j - k + kappa - nu - Max(0,k - nu),2)*((n - q)*psiexp(-j + k,2) + Max(0,k - nu)*(-1 + psiexp(1,2))*psiexp(-j + k,2) + (1 - n + q)*psiexp(1 - j + k,2) + (-1 - j + k - n + q)*psiexp(1 + Max(0,k - nu),2) + (j - k + n - q)*psiexp(2 + Max(0,k - nu),2)) - vEG2G2G2G2cAWBD*b[kappa]*Power2(1/(psiexp(1,1) - psiexp(1,2)))*psiexp(-k + nu,2)*psiexp(-Max(0,k - nu),1)*((-1 - j + k - n + q)*psiexp(1 - j + k,2)*psiexp(1 + Max(0,k - nu),1) + Max(0,k - nu)*(psiexp(1,1) - psiexp(1,2))*psiexp(-j + k,1)*psiexp(1 + Max(0,k - nu),2) + (1 - n + q)*psiexp(1 - j + k,1)*psiexp(1 + Max(0,k - nu),2) + (j - k + n - q)*psiexp(-j + k,2)*psiexp(2 + Max(0,k - nu),1) + (n - q)*psiexp(-j + k,1)*psiexp(2 + Max(0,k - nu),2)) + vEG2G2G2G2cADBH*b[kappa]*Power2(1/(psiexp(1,1) - psiexp(1,2)))*psiexp(-k + nu,1)*psiexp(j - k + kappa - nu - Max(0,k - nu),2)*(Max(0,k - nu)*(psiexp(1,1) - psiexp(1,2))*psiexp(-j + k,2)*psiexp(1 + Max(0,k - nu),1) + (-1 + n - q)*psiexp(1 - j + k,2)*psiexp(1 + Max(0,k - nu),1) + (1 + j - k + n - q)*psiexp(1 - j + k,1)*psiexp(1 + Max(0,k - nu),2) + (-n + q)*psiexp(-j + k,2)*psiexp(2 + Max(0,k - nu),1) + (-j + k - n + q)*psiexp(-j + k,1)*psiexp(2 + Max(0,k - nu),2)) + vEG2G2G2G2cAWAH*b[kappa]*Power2(1/(-1 + psiexp(2,1)))*psiexp(j - k + kappa - nu - 2*Max(0,k - nu),1)*((-n + q)*psiexp(2*(-j + k),1) - Max(0,k - nu)*(-1 + psiexp(2,1))*psiexp(2*(-j + k),1) + (-1 + n - q)*psiexp(2*(1 - j + k),1) + (1 + j - k + n - q)*psiexp(2*(1 + Max(0,k - nu)),1) + (-j + k - n + q)*psiexp(2*(2 + Max(0,k - nu)),1)) - vEG2G2G2G2cAWBDCH*b[kappa]*Power2(1/(psiexp(1,2) - psiexp(1,1)*psiexp(1,3)))*psiexp(-j - k + nu,2)*psiexp(-j - Max(0,k - nu),1)*psiexp(-k + kappa - nu - Max(0,k - nu),3)*(Max(0,k - nu)*(-psiexp(1,2) + psiexp(1,1)*psiexp(1,3))*psiexp(k,1)*psiexp(k,3)*psiexp(1 + j + Max(0,k - nu),2) + (1 - n + q)*psiexp(1 + k,1)*psiexp(1 + k,3)*psiexp(1 + j + Max(0,k - nu),2) + (-1 - j + k - n + q)*psiexp(1 + k,2)*psiexp(1 + j + Max(0,k - nu),1)*psiexp(1 + j + Max(0,k - nu),3) + (n - q)*psiexp(k,1)*psiexp(k,3)*psiexp(2 + j + Max(0,k - nu),2) + (j - k + n - q)*psiexp(k,2)*psiexp(2 + j + Max(0,k - nu),1)*psiexp(2 + j + Max(0,k - nu),3)) - vEG2G2G2G2cAWAHBD*b[kappa]*Power2(1/(-psiexp(1,2) + psiexp(2,1)))*psiexp(-k + nu,2)*psiexp(j - k + kappa - nu - 2*Max(0,k - nu),1)*(Max(0,k - nu)*(-psiexp(1,2) + psiexp(2,1))*psiexp(2*(-j + k),1)*psiexp(1 + Max(0,k - nu),2) + (1 - n + q)*psiexp(2 - 2*j + 2*k,1)*psiexp(1 + Max(0,k - nu),2) + (-1 - j + k - n + q)*psiexp(1 - j + k,2)*psiexp(2*(1 + Max(0,k - nu)),1) + (n - q)*psiexp(2*(-j + k),1)*psiexp(2 + Max(0,k - nu),2) + (j - k + n - q)*psiexp(-j + k,2)*psiexp(4 + 2*Max(0,k - nu),1));
					}
					sumInuP1 = sumIkappa;
				}
				sumInu += sumInuP1*b[nu];
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += MMjkP1/(n - q);
	}
}

for (int j = 1; j <= -4 + q; j++)
{
	for (int k = Max(1,2 + j); k <= -2 + q; k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = 2 + j; nu <= -2 + q; nu++)
			{
				scalar sumInuP1;
				{
					scalar sumIkappa = 0;
					for (int kappa = Max(2 + k,2 + nu); kappa <= Min(-j + k + nu,q); kappa++)
					{
						sumIkappa += -(vEG2G2G2G2c*(-1 + kappa - nu - Max(0,k - nu))*(kappa - 2*n - nu + 2*q + Max(0,k - nu))*b[kappa])/2. - (vEG2G2G2G2cADAH*(-1 + kappa - nu - Max(0,k - nu))*(kappa - 2*n - nu + 2*q + Max(0,k - nu))*b[kappa]*psiexp(-k + kappa,1))/2. - (vEG2G2G2G2cAWAD*(-kappa + 2*n + nu - 2*q - Max(0,k - nu))*(1 - kappa + nu + Max(0,k - nu))*b[kappa]*psiexp(-j + nu,1))/2. + vEG2G2G2G2cAD*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-k + nu,1)*((kappa - n - nu + q)*psiexp(kappa - nu,1) + (1 - kappa + n + nu - q)*psiexp(1 + kappa - nu,1) + (-1 + n - q)*psiexp(1 + Max(0,k - nu),1) + Max(0,k - nu)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,k - nu),1) + (-n + q)*psiexp(2 + Max(0,k - nu),1)) + vEG2G2G2G2cAWADAH*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-j + nu - Max(0,k - nu),1)*((-n + q)*psiexp(kappa - nu,1) - Max(0,k - nu)*(-1 + psiexp(1,1))*psiexp(kappa - nu,1) + (-1 + n - q)*psiexp(1 + kappa - nu,1) + (1 - kappa + n + nu - q)*psiexp(1 + Max(0,k - nu),1) + (kappa - n - nu + q)*psiexp(2 + Max(0,k - nu),1)) + vEG2G2G2G2cAW*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-j + k - kappa + nu - Max(0,k - nu),1)*((-n + q)*psiexp(kappa - nu,1) - Max(0,k - nu)*(-1 + psiexp(1,1))*psiexp(kappa - nu,1) + (-1 + n - q)*psiexp(1 + kappa - nu,1) + (1 - kappa + n + nu - q)*psiexp(1 + Max(0,k - nu),1) + (kappa - n - nu + q)*psiexp(2 + Max(0,k - nu),1)) + vEG2G2G2G2cAWBDBH*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-k + kappa,2)*psiexp(-j + k - kappa + nu - Max(0,k - nu),1)*((-n + q)*psiexp(kappa - nu,1) - Max(0,k - nu)*(-1 + psiexp(1,1))*psiexp(kappa - nu,1) + (-1 + n - q)*psiexp(1 + kappa - nu,1) + (1 - kappa + n + nu - q)*psiexp(1 + Max(0,k - nu),1) + (kappa - n - nu + q)*psiexp(2 + Max(0,k - nu),1)) + vEG2G2G2G2cAH*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-Max(0,k - nu),1)*((-n + q)*psiexp(kappa - nu,1) - Max(0,k - nu)*(-1 + psiexp(1,1))*psiexp(kappa - nu,1) + (-1 + n - q)*psiexp(1 + kappa - nu,1) + (1 - kappa + n + nu - q)*psiexp(1 + Max(0,k - nu),1) + (kappa - n - nu + q)*psiexp(2 + Max(0,k - nu),1)) - vEG2G2G2G2cAWADBH*b[kappa]*Power2(1/(-1 + psiexp(1,2)))*psiexp(-j + nu,1)*psiexp(-Max(0,k - nu),2)*((n - q)*psiexp(kappa - nu,2) + Max(0,k - nu)*(-1 + psiexp(1,2))*psiexp(kappa - nu,2) + (1 - n + q)*psiexp(1 + kappa - nu,2) + (-1 + kappa - n - nu + q)*psiexp(1 + Max(0,k - nu),2) + (-kappa + n + nu - q)*psiexp(2 + Max(0,k - nu),2)) - vEG2G2G2G2cADBH*b[kappa]*Power2(1/(psiexp(1,1) - psiexp(1,2)))*psiexp(-k + nu,1)*psiexp(-Max(0,k - nu),2)*(-(Max(0,k - nu)*(psiexp(1,1) - psiexp(1,2))*psiexp(kappa - nu,2)*psiexp(1 + Max(0,k - nu),1)) + (1 - n + q)*psiexp(1 + kappa - nu,2)*psiexp(1 + Max(0,k - nu),1) + (-1 + kappa - n - nu + q)*psiexp(1 + kappa - nu,1)*psiexp(1 + Max(0,k - nu),2) + (n - q)*psiexp(kappa - nu,2)*psiexp(2 + Max(0,k - nu),1) + (-kappa + n + nu - q)*psiexp(kappa - nu,1)*psiexp(2 + Max(0,k - nu),2)) + vEG2G2G2G2cAWBD*b[kappa]*Power2(1/(psiexp(1,1) - psiexp(1,2)))*psiexp(-k + nu,2)*psiexp(-j + k - kappa + nu - Max(0,k - nu),1)*((1 - kappa + n + nu - q)*psiexp(1 + kappa - nu,2)*psiexp(1 + Max(0,k - nu),1) + Max(0,k - nu)*(-psiexp(1,1) + psiexp(1,2))*psiexp(kappa - nu,1)*psiexp(1 + Max(0,k - nu),2) + (-1 + n - q)*psiexp(1 + kappa - nu,1)*psiexp(1 + Max(0,k - nu),2) + (kappa - n - nu + q)*psiexp(kappa - nu,2)*psiexp(2 + Max(0,k - nu),1) + (-n + q)*psiexp(kappa - nu,1)*psiexp(2 + Max(0,k - nu),2)) + vEG2G2G2G2cAWAH*b[kappa]*Power2(1/(-1 + psiexp(2,1)))*psiexp(-j + k - kappa + nu - 2*Max(0,k - nu),1)*((-n + q)*psiexp(2*(kappa - nu),1) - Max(0,k - nu)*(-1 + psiexp(2,1))*psiexp(2*(kappa - nu),1) + (-1 + n - q)*psiexp(2*(1 + kappa - nu),1) + (1 - kappa + n + nu - q)*psiexp(2*(1 + Max(0,k - nu)),1) + (kappa - n - nu + q)*psiexp(2*(2 + Max(0,k - nu)),1)) + vEG2G2G2G2cAWBDCH*b[kappa]*Power2(1/(psiexp(1,2) - psiexp(1,1)*psiexp(1,3)))*psiexp(-k,2)*psiexp(-j + k - kappa - Max(0,k - nu),1)*psiexp(-nu - Max(0,k - nu),3)*(Max(0,k - nu)*(psiexp(1,2) - psiexp(1,1)*psiexp(1,3))*psiexp(kappa,1)*psiexp(kappa,3)*psiexp(1 + nu + Max(0,k - nu),2) + (-1 + n - q)*psiexp(1 + kappa,1)*psiexp(1 + kappa,3)*psiexp(1 + nu + Max(0,k - nu),2) + (1 - kappa + n + nu - q)*psiexp(1 + kappa,2)*psiexp(1 + nu + Max(0,k - nu),1)*psiexp(1 + nu + Max(0,k - nu),3) + (-n + q)*psiexp(kappa,1)*psiexp(kappa,3)*psiexp(2 + nu + Max(0,k - nu),2) + (kappa - n - nu + q)*psiexp(kappa,2)*psiexp(2 + nu + Max(0,k - nu),1)*psiexp(2 + nu + Max(0,k - nu),3)) + vEG2G2G2G2cAWAHBD*b[kappa]*Power2(1/(-psiexp(1,2) + psiexp(2,1)))*psiexp(-k + nu,2)*psiexp(-j + k - kappa + nu - 2*Max(0,k - nu),1)*(Max(0,k - nu)*(psiexp(1,2) - psiexp(2,1))*psiexp(2*kappa - 2*nu,1)*psiexp(1 + Max(0,k - nu),2) + (-1 + n - q)*psiexp(2*(1 + kappa - nu),1)*psiexp(1 + Max(0,k - nu),2) + (1 - kappa + n + nu - q)*psiexp(1 + kappa - nu,2)*psiexp(2*(1 + Max(0,k - nu)),1) + (-n + q)*psiexp(2*(kappa - nu),1)*psiexp(2 + Max(0,k - nu),2) + (kappa - n - nu + q)*psiexp(kappa - nu,2)*psiexp(4 + 2*Max(0,k - nu),1));
					}
					sumInuP1 = sumIkappa;
				}
				sumInu += sumInuP1*b[nu];
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += MMjkP1/(n - q);
	}
}

for (int j = 1; j <= -3 + q; j++)
{
	for (int k = Max(1,2 + j); k <= -1 + q; k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = 3 + j; nu <= q; nu++)
			{
				scalar sumInuP1;
				{
					scalar sumIkappa = 0;
					for (int kappa = 1 + k; kappa <= Min(-2 - j + k + nu,q); kappa++)
					{
						sumIkappa += -(vEG2G2G2G2c*(j - k + 2*n - 2*q - Max(0,kappa - nu))*(1 + j - k + Max(0,kappa - nu))*b[kappa])/2. - (vEG2G2G2G2cAD*(j - k + 2*n - 2*q - Max(0,kappa - nu))*(1 + j - k + Max(0,kappa - nu))*b[kappa]*psiexp(-k + kappa,1))/2. - (vEG2G2G2G2cAWADAH*(j - k + 2*n - 2*q - Max(0,kappa - nu))*(1 + j - k + Max(0,kappa - nu))*b[kappa]*psiexp(-j + nu,1))/2. - (vEG2G2G2G2cAWAH*(j - k + 2*n - 2*q - Max(0,kappa - nu))*(1 + j - k + Max(0,kappa - nu))*b[kappa]*psiexp(-j + k - kappa + nu,1))/2. - (vEG2G2G2G2cAWAHBD*(j - k + 2*n - 2*q - Max(0,kappa - nu))*(1 + j - k + Max(0,kappa - nu))*b[kappa]*psiexp(-k + kappa,2)*psiexp(-j + k - kappa + nu,1))/2. - vEG2G2G2G2cAWAD*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-k + kappa - Max(0,kappa - nu),1)*((n - q)*psiexp(-j + k,1) + Max(0,kappa - nu)*(-1 + psiexp(1,1))*psiexp(-j + k,1) + (1 - n + q)*psiexp(1 - j + k,1) + (-1 - j + k - n + q)*psiexp(1 + Max(0,kappa - nu),1) + (j - k + n - q)*psiexp(2 + Max(0,kappa - nu),1)) - vEG2G2G2G2cAW*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-Max(0,kappa - nu),1)*((n - q)*psiexp(-j + k,1) + Max(0,kappa - nu)*(-1 + psiexp(1,1))*psiexp(-j + k,1) + (1 - n + q)*psiexp(1 - j + k,1) + (-1 - j + k - n + q)*psiexp(1 + Max(0,kappa - nu),1) + (j - k + n - q)*psiexp(2 + Max(0,kappa - nu),1)) + vEG2G2G2G2cADAH*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-k + nu,1)*((-j + k - n + q)*psiexp(-j + k,1) + (1 + j - k + n - q)*psiexp(1 - j + k,1) + (-1 + n - q)*psiexp(1 + Max(0,kappa - nu),1) + Max(0,kappa - nu)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,kappa - nu),1) + (-n + q)*psiexp(2 + Max(0,kappa - nu),1)) + vEG2G2G2G2cAH*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-kappa + nu,1)*((-j + k - n + q)*psiexp(-j + k,1) + (1 + j - k + n - q)*psiexp(1 - j + k,1) + (-1 + n - q)*psiexp(1 + Max(0,kappa - nu),1) + Max(0,kappa - nu)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,kappa - nu),1) + (-n + q)*psiexp(2 + Max(0,kappa - nu),1)) + vEG2G2G2G2cAWBD*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-k + kappa,2)*psiexp(-Max(0,kappa - nu),1)*((-n + q)*psiexp(-j + k,1) - Max(0,kappa - nu)*(-1 + psiexp(1,1))*psiexp(-j + k,1) + (-1 + n - q)*psiexp(1 - j + k,1) + (1 + j - k + n - q)*psiexp(1 + Max(0,kappa - nu),1) + (-j + k - n + q)*psiexp(2 + Max(0,kappa - nu),1)) + vEG2G2G2G2cADBH*b[kappa]*Power2(1/(-1 + psiexp(1,2)))*psiexp(-k + kappa,1)*psiexp(-kappa + nu,2)*((-j + k - n + q)*psiexp(-j + k,2) + (1 + j - k + n - q)*psiexp(1 - j + k,2) + (-1 + n - q)*psiexp(1 + Max(0,kappa - nu),2) + Max(0,kappa - nu)*(-1 + psiexp(1,2))*psiexp(1 + Max(0,kappa - nu),2) + (-n + q)*psiexp(2 + Max(0,kappa - nu),2)) - vEG2G2G2G2cAWBDBH*b[kappa]*Power2(1/(psiexp(1,1) - psiexp(1,2)))*psiexp(-k + nu,2)*psiexp(-Max(0,kappa - nu),1)*((-1 - j + k - n + q)*psiexp(1 - j + k,2)*psiexp(1 + Max(0,kappa - nu),1) + Max(0,kappa - nu)*(psiexp(1,1) - psiexp(1,2))*psiexp(-j + k,1)*psiexp(1 + Max(0,kappa - nu),2) + (1 - n + q)*psiexp(1 - j + k,1)*psiexp(1 + Max(0,kappa - nu),2) + (j - k + n - q)*psiexp(-j + k,2)*psiexp(2 + Max(0,kappa - nu),1) + (n - q)*psiexp(-j + k,1)*psiexp(2 + Max(0,kappa - nu),2)) + vEG2G2G2G2cAWADBH*b[kappa]*Power2(1/(psiexp(1,1) - psiexp(1,2)))*psiexp(-kappa + nu,2)*psiexp(-k + kappa - Max(0,kappa - nu),1)*((1 + j - k + n - q)*psiexp(1 - j + k,2)*psiexp(1 + Max(0,kappa - nu),1) + Max(0,kappa - nu)*(-psiexp(1,1) + psiexp(1,2))*psiexp(-j + k,1)*psiexp(1 + Max(0,kappa - nu),2) + (-1 + n - q)*psiexp(1 - j + k,1)*psiexp(1 + Max(0,kappa - nu),2) + (-j + k - n + q)*psiexp(-j + k,2)*psiexp(2 + Max(0,kappa - nu),1) + (-n + q)*psiexp(-j + k,1)*psiexp(2 + Max(0,kappa - nu),2)) - vEG2G2G2G2cAWBDCH*b[kappa]*Power2(1/(psiexp(1,1) - psiexp(1,3)))*psiexp(-k + kappa,2)*psiexp(-kappa + nu,3)*psiexp(-Max(0,kappa - nu),1)*((-1 - j + k - n + q)*psiexp(1 - j + k,3)*psiexp(1 + Max(0,kappa - nu),1) + Max(0,kappa - nu)*(psiexp(1,1) - psiexp(1,3))*psiexp(-j + k,1)*psiexp(1 + Max(0,kappa - nu),3) + (1 - n + q)*psiexp(1 - j + k,1)*psiexp(1 + Max(0,kappa - nu),3) + (j - k + n - q)*psiexp(-j + k,3)*psiexp(2 + Max(0,kappa - nu),1) + (n - q)*psiexp(-j + k,1)*psiexp(2 + Max(0,kappa - nu),3));
					}
					sumInuP1 = sumIkappa;
				}
				sumInu += sumInuP1*b[nu];
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += MMjkP1/(n - q);
	}
}

for (int j = 1; j <= -3 + q; j++)
{
	for (int k = Max(1,3 + j); k <= q; k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = 3 + j; nu <= q; nu++)
			{
				scalar sumInuP1;
				{
					scalar sumIkappa = 0;
					for (int kappa = Max3(0,2 + j,2 + j + k - nu); kappa <= -1 + k; kappa++)
					{
						sumIkappa += -(vEG2G2G2G2c*(j - kappa + 2*n - 2*q - Max(0,k - nu))*(1 + j - kappa + Max(0,k - nu))*b[kappa])/2. - (vEG2G2G2G2cAD*(-1 - j + kappa - Max(0,k - nu))*(-j + kappa - 2*n + 2*q + Max(0,k - nu))*b[kappa]*psiexp(k - kappa,1))/2. - (vEG2G2G2G2cAWADAH*(j - kappa + 2*n - 2*q - Max(0,k - nu))*(1 + j - kappa + Max(0,k - nu))*b[kappa]*psiexp(-j + nu,1))/2. - (vEG2G2G2G2cAWAH*(j - kappa + 2*n - 2*q - Max(0,k - nu))*(1 + j - kappa + Max(0,k - nu))*b[kappa]*psiexp(-j - k + kappa + nu,1))/2. - (vEG2G2G2G2cAWAHBD*(j - kappa + 2*n - 2*q - Max(0,k - nu))*(1 + j - kappa + Max(0,k - nu))*b[kappa]*psiexp(k - kappa,2)*psiexp(-j - k + kappa + nu,1))/2. - vEG2G2G2G2cAWAD*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(k - kappa - Max(0,k - nu),1)*((n - q)*psiexp(-j + kappa,1) + Max(0,k - nu)*(-1 + psiexp(1,1))*psiexp(-j + kappa,1) + (1 - n + q)*psiexp(1 - j + kappa,1) + (-1 - j + kappa - n + q)*psiexp(1 + Max(0,k - nu),1) + (j - kappa + n - q)*psiexp(2 + Max(0,k - nu),1)) - vEG2G2G2G2cAW*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-Max(0,k - nu),1)*((n - q)*psiexp(-j + kappa,1) + Max(0,k - nu)*(-1 + psiexp(1,1))*psiexp(-j + kappa,1) + (1 - n + q)*psiexp(1 - j + kappa,1) + (-1 - j + kappa - n + q)*psiexp(1 + Max(0,k - nu),1) + (j - kappa + n - q)*psiexp(2 + Max(0,k - nu),1)) + vEG2G2G2G2cAH*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-k + nu,1)*((-j + kappa - n + q)*psiexp(-j + kappa,1) + (1 + j - kappa + n - q)*psiexp(1 - j + kappa,1) + (-1 + n - q)*psiexp(1 + Max(0,k - nu),1) + Max(0,k - nu)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,k - nu),1) + (-n + q)*psiexp(2 + Max(0,k - nu),1)) + vEG2G2G2G2cADAH*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-kappa + nu,1)*((-j + kappa - n + q)*psiexp(-j + kappa,1) + (1 + j - kappa + n - q)*psiexp(1 - j + kappa,1) + (-1 + n - q)*psiexp(1 + Max(0,k - nu),1) + Max(0,k - nu)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,k - nu),1) + (-n + q)*psiexp(2 + Max(0,k - nu),1)) + vEG2G2G2G2cAWBD*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(k - kappa,2)*psiexp(-Max(0,k - nu),1)*((-n + q)*psiexp(-j + kappa,1) - Max(0,k - nu)*(-1 + psiexp(1,1))*psiexp(-j + kappa,1) + (-1 + n - q)*psiexp(1 - j + kappa,1) + (1 + j - kappa + n - q)*psiexp(1 + Max(0,k - nu),1) + (-j + kappa - n + q)*psiexp(2 + Max(0,k - nu),1)) + vEG2G2G2G2cADBH*b[kappa]*Power2(1/(-1 + psiexp(1,2)))*psiexp(k - kappa,1)*psiexp(-k + nu,2)*((-j + kappa - n + q)*psiexp(-j + kappa,2) + (1 + j - kappa + n - q)*psiexp(1 - j + kappa,2) + (-1 + n - q)*psiexp(1 + Max(0,k - nu),2) + Max(0,k - nu)*(-1 + psiexp(1,2))*psiexp(1 + Max(0,k - nu),2) + (-n + q)*psiexp(2 + Max(0,k - nu),2)) - vEG2G2G2G2cAWBDBH*b[kappa]*Power2(1/(psiexp(1,1) - psiexp(1,2)))*psiexp(-kappa + nu,2)*psiexp(-Max(0,k - nu),1)*((-1 - j + kappa - n + q)*psiexp(1 - j + kappa,2)*psiexp(1 + Max(0,k - nu),1) + Max(0,k - nu)*(psiexp(1,1) - psiexp(1,2))*psiexp(-j + kappa,1)*psiexp(1 + Max(0,k - nu),2) + (1 - n + q)*psiexp(1 - j + kappa,1)*psiexp(1 + Max(0,k - nu),2) + (j - kappa + n - q)*psiexp(-j + kappa,2)*psiexp(2 + Max(0,k - nu),1) + (n - q)*psiexp(-j + kappa,1)*psiexp(2 + Max(0,k - nu),2)) + vEG2G2G2G2cAWADBH*b[kappa]*Power2(1/(psiexp(1,1) - psiexp(1,2)))*psiexp(-k + nu,2)*psiexp(k - kappa - Max(0,k - nu),1)*((1 + j - kappa + n - q)*psiexp(1 - j + kappa,2)*psiexp(1 + Max(0,k - nu),1) + Max(0,k - nu)*(-psiexp(1,1) + psiexp(1,2))*psiexp(-j + kappa,1)*psiexp(1 + Max(0,k - nu),2) + (-1 + n - q)*psiexp(1 - j + kappa,1)*psiexp(1 + Max(0,k - nu),2) + (-j + kappa - n + q)*psiexp(-j + kappa,2)*psiexp(2 + Max(0,k - nu),1) + (-n + q)*psiexp(-j + kappa,1)*psiexp(2 + Max(0,k - nu),2)) - vEG2G2G2G2cAWBDCH*b[kappa]*Power2(1/(psiexp(1,1) - psiexp(1,3)))*psiexp(k - kappa,2)*psiexp(-k + nu,3)*psiexp(-Max(0,k - nu),1)*((-1 - j + kappa - n + q)*psiexp(1 - j + kappa,3)*psiexp(1 + Max(0,k - nu),1) + Max(0,k - nu)*(psiexp(1,1) - psiexp(1,3))*psiexp(-j + kappa,1)*psiexp(1 + Max(0,k - nu),3) + (1 - n + q)*psiexp(1 - j + kappa,1)*psiexp(1 + Max(0,k - nu),3) + (j - kappa + n - q)*psiexp(-j + kappa,3)*psiexp(2 + Max(0,k - nu),1) + (n - q)*psiexp(-j + kappa,1)*psiexp(2 + Max(0,k - nu),3));
					}
					sumInuP1 = sumIkappa;
				}
				sumInu += sumInuP1*b[nu];
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += MMjkP1/(n - q);
	}
}

for (int j = 1; j <= -2 + q; j++)
{
	for (int k = Max(2,2 + j); k <= q; k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = 2 + j; nu <= Min(j + k,q); nu++)
			{
				scalar sumInuP1;
				{
					scalar sumIkappa = 0;
					for (int kappa = 0; kappa <= j + k - nu; kappa++)
					{
						sumIkappa += -(vEG2G2G2G2c*(j - k + 2*n - 2*q - Max(0,k - nu))*(1 + j - k + Max(0,k - nu))*b[kappa])/2. - (vEG2G2G2G2cAWAD*(-1 - j + k - Max(0,k - nu))*(-j + k - 2*n + 2*q + Max(0,k - nu))*b[kappa]*psiexp(k - kappa,1))/2. - (vEG2G2G2G2cADAH*(j - k + 2*n - 2*q - Max(0,k - nu))*(1 + j - k + Max(0,k - nu))*b[kappa]*psiexp(-j + nu,1))/2. - vEG2G2G2G2cAD*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-Max(0,k - nu),1)*((n - q)*psiexp(-j + k,1) + Max(0,k - nu)*(-1 + psiexp(1,1))*psiexp(-j + k,1) + (1 - n + q)*psiexp(1 - j + k,1) + (-1 - j + k - n + q)*psiexp(1 + Max(0,k - nu),1) + (j - k + n - q)*psiexp(2 + Max(0,k - nu),1)) + vEG2G2G2G2cAW*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(j - kappa,1)*((-j + k - n + q)*psiexp(-j + k,1) + (1 + j - k + n - q)*psiexp(1 - j + k,1) + (-1 + n - q)*psiexp(1 + Max(0,k - nu),1) + Max(0,k - nu)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,k - nu),1) + (-n + q)*psiexp(2 + Max(0,k - nu),1)) + vEG2G2G2G2cAWBDBH*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(j - kappa,1)*psiexp(-j + nu,2)*((-j + k - n + q)*psiexp(-j + k,1) + (1 + j - k + n - q)*psiexp(1 - j + k,1) + (-1 + n - q)*psiexp(1 + Max(0,k - nu),1) + Max(0,k - nu)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,k - nu),1) + (-n + q)*psiexp(2 + Max(0,k - nu),1)) + vEG2G2G2G2cAH*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-k + nu,1)*((-j + k - n + q)*psiexp(-j + k,1) + (1 + j - k + n - q)*psiexp(1 - j + k,1) + (-1 + n - q)*psiexp(1 + Max(0,k - nu),1) + Max(0,k - nu)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,k - nu),1) + (-n + q)*psiexp(2 + Max(0,k - nu),1)) + vEG2G2G2G2cAWADAH*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-kappa + nu,1)*((-j + k - n + q)*psiexp(-j + k,1) + (1 + j - k + n - q)*psiexp(1 - j + k,1) + (-1 + n - q)*psiexp(1 + Max(0,k - nu),1) + Max(0,k - nu)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,k - nu),1) + (-n + q)*psiexp(2 + Max(0,k - nu),1)) + vEG2G2G2G2cAWADBH*b[kappa]*Power2(1/(-1 + psiexp(1,2)))*psiexp(k - kappa,1)*psiexp(-k + nu,2)*((-j + k - n + q)*psiexp(-j + k,2) + (1 + j - k + n - q)*psiexp(1 - j + k,2) + (-1 + n - q)*psiexp(1 + Max(0,k - nu),2) + Max(0,k - nu)*(-1 + psiexp(1,2))*psiexp(1 + Max(0,k - nu),2) + (-n + q)*psiexp(2 + Max(0,k - nu),2)) - vEG2G2G2G2cADBH*b[kappa]*Power2(1/(psiexp(1,1) - psiexp(1,2)))*psiexp(-k + nu,2)*psiexp(-Max(0,k - nu),1)*((-1 - j + k - n + q)*psiexp(1 - j + k,2)*psiexp(1 + Max(0,k - nu),1) + Max(0,k - nu)*(psiexp(1,1) - psiexp(1,2))*psiexp(-j + k,1)*psiexp(1 + Max(0,k - nu),2) + (1 - n + q)*psiexp(1 - j + k,1)*psiexp(1 + Max(0,k - nu),2) + (j - k + n - q)*psiexp(-j + k,2)*psiexp(2 + Max(0,k - nu),1) + (n - q)*psiexp(-j + k,1)*psiexp(2 + Max(0,k - nu),2)) - vEG2G2G2G2cAWBD*b[kappa]*Power2(1/(psiexp(1,1) - psiexp(1,2)))*psiexp(j - kappa,1)*psiexp(-Max(0,k - nu),2)*(-(Max(0,k - nu)*(psiexp(1,1) - psiexp(1,2))*psiexp(-j + k,2)*psiexp(1 + Max(0,k - nu),1)) + (1 - n + q)*psiexp(1 - j + k,2)*psiexp(1 + Max(0,k - nu),1) + (-1 - j + k - n + q)*psiexp(1 - j + k,1)*psiexp(1 + Max(0,k - nu),2) + (n - q)*psiexp(-j + k,2)*psiexp(2 + Max(0,k - nu),1) + (j - k + n - q)*psiexp(-j + k,1)*psiexp(2 + Max(0,k - nu),2)) + vEG2G2G2G2cAWAH*b[kappa]*Power2(1/(-1 + psiexp(2,1)))*psiexp(j - k - kappa + nu,1)*((-j + k - n + q)*psiexp(2*(-j + k),1) + (1 + j - k + n - q)*psiexp(2*(1 - j + k),1) + (-1 + n - q)*psiexp(2*(1 + Max(0,k - nu)),1) + Max(0,k - nu)*(-1 + psiexp(2,1))*psiexp(2*(1 + Max(0,k - nu)),1) + (-n + q)*psiexp(2*(2 + Max(0,k - nu)),1)) + vEG2G2G2G2cAWBDCH*b[kappa]*Power2(1/(psiexp(1,2) - psiexp(1,1)*psiexp(1,3)))*psiexp(-kappa,1)*psiexp(-j - k + nu,3)*psiexp(-j - Max(0,k - nu),2)*((1 + j - k + n - q)*psiexp(1 + k,1)*psiexp(1 + k,3)*psiexp(1 + j + Max(0,k - nu),2) + Max(0,k - nu)*(-(psiexp(1,2)*psiexp(1 + j,1)) + psiexp(1,3)*psiexp(2 + j,1))*psiexp(k,2)*psiexp(Max(0,k - nu),1)*psiexp(1 + j + Max(0,k - nu),3) + (-1 + n - q)*psiexp(1 + k,2)*psiexp(1 + j + Max(0,k - nu),1)*psiexp(1 + j + Max(0,k - nu),3) + (-j + k - n + q)*psiexp(k,1)*psiexp(k,3)*psiexp(2 + j + Max(0,k - nu),2) + (-n + q)*psiexp(k,2)*psiexp(2 + j + Max(0,k - nu),1)*psiexp(2 + j + Max(0,k - nu),3)) - vEG2G2G2G2cAWAHBD*b[kappa]*Power2(1/(-psiexp(1,2) + psiexp(2,1)))*psiexp(j - k - kappa + nu,1)*psiexp(-Max(0,k - nu),2)*((-1 - j + k - n + q)*psiexp(2 - 2*j + 2*k,1)*psiexp(1 + Max(0,k - nu),2) - Max(0,k - nu)*(-psiexp(1,2) + psiexp(2,1))*psiexp(-j + k,2)*psiexp(2*(1 + Max(0,k - nu)),1) + (1 - n + q)*psiexp(1 - j + k,2)*psiexp(2*(1 + Max(0,k - nu)),1) + (j - k + n - q)*psiexp(2*(-j + k),1)*psiexp(2 + Max(0,k - nu),2) + (n - q)*psiexp(-j + k,2)*psiexp(4 + 2*Max(0,k - nu),1));
					}
					sumInuP1 = sumIkappa;
				}
				sumInu += sumInuP1*b[nu];
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += MMjkP1/(n - q);
	}
}

for (int j = 1; j <= -3 + q; j++)
{
	for (int k = Max(2,2 + j); k <= q; k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = 3 + j; nu <= q; nu++)
			{
				scalar sumInuP1;
				{
					scalar sumIkappa = 0;
					for (int kappa = Max(0,1 + j + k - nu); kappa <= -2 + k; kappa++)
					{
						sumIkappa += -(vEG2G2G2G2c*(j - k + 2*n - 2*q - Max(0,-j + kappa))*(1 + j - k + Max(0,-j + kappa))*b[kappa])/2. - (vEG2G2G2G2cAWAD*(-1 - j + k - Max(0,-j + kappa))*(-j + k - 2*n + 2*q + Max(0,-j + kappa))*b[kappa]*psiexp(k - kappa,1))/2. - (vEG2G2G2G2cADAH*(j - k + 2*n - 2*q - Max(0,-j + kappa))*(1 + j - k + Max(0,-j + kappa))*b[kappa]*psiexp(-j + nu,1))/2. - vEG2G2G2G2cAD*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-Max(0,-j + kappa),1)*((n - q)*psiexp(-j + k,1) + Max(0,-j + kappa)*(-1 + psiexp(1,1))*psiexp(-j + k,1) + (1 - n + q)*psiexp(1 - j + k,1) + (-1 - j + k - n + q)*psiexp(1 + Max(0,-j + kappa),1) + (j - k + n - q)*psiexp(2 + Max(0,-j + kappa),1)) + vEG2G2G2G2cAW*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(j - kappa,1)*((-j + k - n + q)*psiexp(-j + k,1) + (1 + j - k + n - q)*psiexp(1 - j + k,1) + (-1 + n - q)*psiexp(1 + Max(0,-j + kappa),1) + Max(0,-j + kappa)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,-j + kappa),1) + (-n + q)*psiexp(2 + Max(0,-j + kappa),1)) + vEG2G2G2G2cAWBDBH*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(j - kappa,1)*psiexp(-j + nu,2)*((-j + k - n + q)*psiexp(-j + k,1) + (1 + j - k + n - q)*psiexp(1 - j + k,1) + (-1 + n - q)*psiexp(1 + Max(0,-j + kappa),1) + Max(0,-j + kappa)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,-j + kappa),1) + (-n + q)*psiexp(2 + Max(0,-j + kappa),1)) + vEG2G2G2G2cAH*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-k + nu,1)*((-j + k - n + q)*psiexp(-j + k,1) + (1 + j - k + n - q)*psiexp(1 - j + k,1) + (-1 + n - q)*psiexp(1 + Max(0,-j + kappa),1) + Max(0,-j + kappa)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,-j + kappa),1) + (-n + q)*psiexp(2 + Max(0,-j + kappa),1)) + vEG2G2G2G2cAWADAH*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-kappa + nu,1)*((-j + k - n + q)*psiexp(-j + k,1) + (1 + j - k + n - q)*psiexp(1 - j + k,1) + (-1 + n - q)*psiexp(1 + Max(0,-j + kappa),1) + Max(0,-j + kappa)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,-j + kappa),1) + (-n + q)*psiexp(2 + Max(0,-j + kappa),1)) + vEG2G2G2G2cAWADBH*b[kappa]*Power2(1/(-1 + psiexp(1,2)))*psiexp(k - kappa,1)*psiexp(-k + nu,2)*((-j + k - n + q)*psiexp(-j + k,2) + (1 + j - k + n - q)*psiexp(1 - j + k,2) + (-1 + n - q)*psiexp(1 + Max(0,-j + kappa),2) + Max(0,-j + kappa)*(-1 + psiexp(1,2))*psiexp(1 + Max(0,-j + kappa),2) + (-n + q)*psiexp(2 + Max(0,-j + kappa),2)) - vEG2G2G2G2cADBH*b[kappa]*Power2(1/(psiexp(1,1) - psiexp(1,2)))*psiexp(-k + nu,2)*psiexp(-Max(0,-j + kappa),1)*((-1 - j + k - n + q)*psiexp(1 - j + k,2)*psiexp(1 + Max(0,-j + kappa),1) + Max(0,-j + kappa)*(psiexp(1,1) - psiexp(1,2))*psiexp(-j + k,1)*psiexp(1 + Max(0,-j + kappa),2) + (1 - n + q)*psiexp(1 - j + k,1)*psiexp(1 + Max(0,-j + kappa),2) + (j - k + n - q)*psiexp(-j + k,2)*psiexp(2 + Max(0,-j + kappa),1) + (n - q)*psiexp(-j + k,1)*psiexp(2 + Max(0,-j + kappa),2)) - vEG2G2G2G2cAWBD*b[kappa]*Power2(1/(psiexp(1,1) - psiexp(1,2)))*psiexp(j - kappa,1)*psiexp(-Max(0,-j + kappa),2)*(-(Max(0,-j + kappa)*(psiexp(1,1) - psiexp(1,2))*psiexp(-j + k,2)*psiexp(1 + Max(0,-j + kappa),1)) + (1 - n + q)*psiexp(1 - j + k,2)*psiexp(1 + Max(0,-j + kappa),1) + (-1 - j + k - n + q)*psiexp(1 - j + k,1)*psiexp(1 + Max(0,-j + kappa),2) + (n - q)*psiexp(-j + k,2)*psiexp(2 + Max(0,-j + kappa),1) + (j - k + n - q)*psiexp(-j + k,1)*psiexp(2 + Max(0,-j + kappa),2)) + vEG2G2G2G2cAWAH*b[kappa]*Power2(1/(-1 + psiexp(2,1)))*psiexp(j - k - kappa + nu,1)*((-j + k - n + q)*psiexp(2*(-j + k),1) + (1 + j - k + n - q)*psiexp(2*(1 - j + k),1) + (-1 + n - q)*psiexp(2*(1 + Max(0,-j + kappa)),1) + Max(0,-j + kappa)*(-1 + psiexp(2,1))*psiexp(2*(1 + Max(0,-j + kappa)),1) + (-n + q)*psiexp(2*(2 + Max(0,-j + kappa)),1)) + vEG2G2G2G2cAWBDCH*b[kappa]*Power2(1/(psiexp(1,2) - psiexp(1,1)*psiexp(1,3)))*psiexp(-kappa,1)*psiexp(-j - k + nu,3)*psiexp(-j - Max(0,-j + kappa),2)*((1 + j - k + n - q)*psiexp(1 + k,1)*psiexp(1 + k,3)*psiexp(1 + j + Max(0,-j + kappa),2) + Max(0,-j + kappa)*(-(psiexp(1,2)*psiexp(1 + j,1)) + psiexp(1,3)*psiexp(2 + j,1))*psiexp(k,2)*psiexp(Max(0,-j + kappa),1)*psiexp(1 + j + Max(0,-j + kappa),3) + (-1 + n - q)*psiexp(1 + k,2)*psiexp(1 + j + Max(0,-j + kappa),1)*psiexp(1 + j + Max(0,-j + kappa),3) + (-j + k - n + q)*psiexp(k,1)*psiexp(k,3)*psiexp(2 + j + Max(0,-j + kappa),2) + (-n + q)*psiexp(k,2)*psiexp(2 + j + Max(0,-j + kappa),1)*psiexp(2 + j + Max(0,-j + kappa),3)) - vEG2G2G2G2cAWAHBD*b[kappa]*Power2(1/(-psiexp(1,2) + psiexp(2,1)))*psiexp(j - k - kappa + nu,1)*psiexp(-Max(0,-j + kappa),2)*((-1 - j + k - n + q)*psiexp(2 - 2*j + 2*k,1)*psiexp(1 + Max(0,-j + kappa),2) - Max(0,-j + kappa)*(-psiexp(1,2) + psiexp(2,1))*psiexp(-j + k,2)*psiexp(2*(1 + Max(0,-j + kappa)),1) + (1 - n + q)*psiexp(1 - j + k,2)*psiexp(2*(1 + Max(0,-j + kappa)),1) + (j - k + n - q)*psiexp(2*(-j + k),1)*psiexp(2 + Max(0,-j + kappa),2) + (n - q)*psiexp(-j + k,2)*psiexp(4 + 2*Max(0,-j + kappa),1));
					}
					sumInuP1 = sumIkappa;
				}
				sumInu += sumInuP1*b[nu];
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += MMjkP1/(n - q);
	}
}

for (int j = 1; j <= -3 + q; j++)
{
	for (int k = Max(1,2 + j); k <= -1 + q; k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = 2 + j; nu <= -1 + q; nu++)
			{
				scalar sumInuP1;
				{
					scalar sumIkappa = 0;
					for (int kappa = Max(1 + k,1 + nu); kappa <= Min(-1 - j + k + nu,q); kappa++)
					{
						sumIkappa += (-kappa + n + nu - q)*vEG2G2G4c*b[kappa] + (-kappa + n + nu - q)*vEG2G2G4cAD*b[kappa]*psiexp(-k + kappa,1) + (-kappa + n + nu - q)*vEG2G2G4cAWAD*b[kappa]*psiexp(-j + nu,1) - (kappa - n - nu + q)*vEG2G2G4cAW*b[kappa]*psiexp(-j + k - kappa + nu,1) + (-kappa + n + nu - q)*vEG2G2G4cAWBD*b[kappa]*psiexp(-k + kappa,2)*psiexp(-j + k - kappa + nu,1);
					}
					sumInuP1 = sumIkappa;
				}
				sumInu += sumInuP1*b[nu];
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += MMjkP1/(n - q);
	}
}

for (int j = 1; j <= -2 + q; j++)
{
	for (int k = Max(1,2 + j); k <= q; k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = 2 + j; nu <= q; nu++)
			{
				sumInu += -(vEG2G4G2c*(j - k + 2*n - 2*q - Max(0,k - nu))*(1 + j - k + Max(0,k - nu))*b[nu])/2. - (vEG2G4G2cAWAH*(j - k + 2*n - 2*q - Max(0,k - nu))*(1 + j - k + Max(0,k - nu))*b[nu]*psiexp(-j + nu,1))/2. - vEG2G4G2cAW*b[nu]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-Max(0,k - nu),1)*((n - q)*psiexp(-j + k,1) + Max(0,k - nu)*(-1 + psiexp(1,1))*psiexp(-j + k,1) + (1 - n + q)*psiexp(1 - j + k,1) + (-1 - j + k - n + q)*psiexp(1 + Max(0,k - nu),1) + (j - k + n - q)*psiexp(2 + Max(0,k - nu),1)) + vEG2G4G2cAH*b[nu]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-k + nu,1)*((-j + k - n + q)*psiexp(-j + k,1) + (1 + j - k + n - q)*psiexp(1 - j + k,1) + (-1 + n - q)*psiexp(1 + Max(0,k - nu),1) + Max(0,k - nu)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,k - nu),1) + (-n + q)*psiexp(2 + Max(0,k - nu),1)) + vEG2G4G2cBH*b[nu]*Power2(1/(-1 + psiexp(1,2)))*psiexp(-k + nu,2)*((-j + k - n + q)*psiexp(-j + k,2) + (1 + j - k + n - q)*psiexp(1 - j + k,2) + (-1 + n - q)*psiexp(1 + Max(0,k - nu),2) + Max(0,k - nu)*(-1 + psiexp(1,2))*psiexp(1 + Max(0,k - nu),2) + (-n + q)*psiexp(2 + Max(0,k - nu),2)) - vEG2G4G2cAWBH*b[nu]*Power2(1/(psiexp(1,1) - psiexp(1,2)))*psiexp(-k + nu,2)*psiexp(-Max(0,k - nu),1)*((-1 - j + k - n + q)*psiexp(1 - j + k,2)*psiexp(1 + Max(0,k - nu),1) + Max(0,k - nu)*(psiexp(1,1) - psiexp(1,2))*psiexp(-j + k,1)*psiexp(1 + Max(0,k - nu),2) + (1 - n + q)*psiexp(1 - j + k,1)*psiexp(1 + Max(0,k - nu),2) + (j - k + n - q)*psiexp(-j + k,2)*psiexp(2 + Max(0,k - nu),1) + (n - q)*psiexp(-j + k,1)*psiexp(2 + Max(0,k - nu),2)) - vEG2G4G2cAWCH*b[nu]*Power2(1/(psiexp(1,1) - psiexp(1,3)))*psiexp(-k + nu,3)*psiexp(-Max(0,k - nu),1)*((-1 - j + k - n + q)*psiexp(1 - j + k,3)*psiexp(1 + Max(0,k - nu),1) + Max(0,k - nu)*(psiexp(1,1) - psiexp(1,3))*psiexp(-j + k,1)*psiexp(1 + Max(0,k - nu),3) + (1 - n + q)*psiexp(1 - j + k,1)*psiexp(1 + Max(0,k - nu),3) + (j - k + n - q)*psiexp(-j + k,3)*psiexp(2 + Max(0,k - nu),1) + (n - q)*psiexp(-j + k,1)*psiexp(2 + Max(0,k - nu),3));
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += (MMjkP1*b[k])/(n - q);
	}
}

for (int j = 1; j <= -2 + q; j++)
{
	for (int k = Max(1,2 + j); k <= q; k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = 2 + j; nu <= q; nu++)
			{
				scalar sumInuP1;
				{
					scalar sumIkappa = 0;
					for (int kappa = Max3(0,1 + j,1 + j + k - nu); kappa <= -1 + k; kappa++)
					{
						sumIkappa += (j - kappa + n - q)*vEG4G2G2c*b[kappa] + (j - kappa + n - q)*vEG4G2G2cAD*b[kappa]*psiexp(k - kappa,1) + (j - kappa + n - q)*vEG4G2G2cBD*b[kappa]*psiexp(k - kappa,2) + (j - kappa + n - q)*vEG4G2G2cADAH*b[kappa]*psiexp(-j + nu,1) + (j - kappa + n - q)*vEG4G2G2cBDBH*b[kappa]*psiexp(-j + nu,2) + (j - kappa + n - q)*vEG4G2G2cAH*b[kappa]*psiexp(-j - k + kappa + nu,1) + (j - kappa + n - q)*vEG4G2G2cAHBD*b[kappa]*psiexp(k - kappa,2)*psiexp(-j - k + kappa + nu,1) + (j - kappa + n - q)*vEG4G2G2cADBH*b[kappa]*psiexp(k - kappa,1)*psiexp(-j - k + kappa + nu,2) + (j - kappa + n - q)*vEG4G2G2cBDCH*b[kappa]*psiexp(k - kappa,2)*psiexp(-j - k + kappa + nu,3);
					}
					sumInuP1 = sumIkappa;
				}
				sumInu += sumInuP1*b[nu];
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += MMjkP1/(n - q);
	}
}

for (int j = 1; j <= -3 + q; j++)
{
	for (int k = 1; k <= -3 + q; k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = 1 + j; nu <= Min(-2 + q,-2 + j - k + q); nu++)
			{
				scalar sumInuP1;
				{
					scalar sumIkappa = 0;
					for (int kappa = Max(2 + nu,2 - j + k + nu); kappa <= q; kappa++)
					{
						sumIkappa += -(vEG2G2G2G2c*(-1 + kappa - nu - Max(0,-j + k))*(kappa - 2*n - nu + 2*q + Max(0,-j + k))*b[kappa])/2. - (vEG2G2G2G2cAWADAH*(-1 + kappa - nu - Max(0,-j + k))*(kappa - 2*n - nu + 2*q + Max(0,-j + k))*b[kappa]*psiexp(-k + kappa,1))/2. - (vEG2G2G2G2cAWAH*(-1 + kappa - nu - Max(0,-j + k))*(kappa - 2*n - nu + 2*q + Max(0,-j + k))*b[kappa]*psiexp(j - k + kappa - nu,1))/2. - (vEG2G2G2G2cAD*(-kappa + 2*n + nu - 2*q - Max(0,-j + k))*(1 - kappa + nu + Max(0,-j + k))*b[kappa]*psiexp(-j + nu,1))/2. - (vEG2G2G2G2cAWAHBD*(-1 + kappa - nu - Max(0,-j + k))*(kappa - 2*n - nu + 2*q + Max(0,-j + k))*b[kappa]*psiexp(j - k + kappa - nu,1)*psiexp(-j + nu,2))/2. + vEG2G2G2G2cAW*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(j - k,1)*((kappa - n - nu + q)*psiexp(kappa - nu,1) + (1 - kappa + n + nu - q)*psiexp(1 + kappa - nu,1) + (-1 + n - q)*psiexp(1 + Max(0,-j + k),1) + Max(0,-j + k)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,-j + k),1) + (-n + q)*psiexp(2 + Max(0,-j + k),1)) + vEG2G2G2G2cAWBD*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(j - k,1)*psiexp(-j + nu,2)*((kappa - n - nu + q)*psiexp(kappa - nu,1) + (1 - kappa + n + nu - q)*psiexp(1 + kappa - nu,1) + (-1 + n - q)*psiexp(1 + Max(0,-j + k),1) + Max(0,-j + k)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,-j + k),1) + (-n + q)*psiexp(2 + Max(0,-j + k),1)) + vEG2G2G2G2cAWAD*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-k + nu,1)*((kappa - n - nu + q)*psiexp(kappa - nu,1) + (1 - kappa + n + nu - q)*psiexp(1 + kappa - nu,1) + (-1 + n - q)*psiexp(1 + Max(0,-j + k),1) + Max(0,-j + k)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,-j + k),1) + (-n + q)*psiexp(2 + Max(0,-j + k),1)) + vEG2G2G2G2cADAH*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-j + nu - Max(0,-j + k),1)*((-n + q)*psiexp(kappa - nu,1) - Max(0,-j + k)*(-1 + psiexp(1,1))*psiexp(kappa - nu,1) + (-1 + n - q)*psiexp(1 + kappa - nu,1) + (1 - kappa + n + nu - q)*psiexp(1 + Max(0,-j + k),1) + (kappa - n - nu + q)*psiexp(2 + Max(0,-j + k),1)) + vEG2G2G2G2cAH*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-Max(0,-j + k),1)*((-n + q)*psiexp(kappa - nu,1) - Max(0,-j + k)*(-1 + psiexp(1,1))*psiexp(kappa - nu,1) + (-1 + n - q)*psiexp(1 + kappa - nu,1) + (1 - kappa + n + nu - q)*psiexp(1 + Max(0,-j + k),1) + (kappa - n - nu + q)*psiexp(2 + Max(0,-j + k),1)) - vEG2G2G2G2cADBH*b[kappa]*Power2(1/(-1 + psiexp(1,2)))*psiexp(-j + nu,1)*psiexp(-Max(0,-j + k),2)*((n - q)*psiexp(kappa - nu,2) + Max(0,-j + k)*(-1 + psiexp(1,2))*psiexp(kappa - nu,2) + (1 - n + q)*psiexp(1 + kappa - nu,2) + (-1 + kappa - n - nu + q)*psiexp(1 + Max(0,-j + k),2) + (-kappa + n + nu - q)*psiexp(2 + Max(0,-j + k),2)) - vEG2G2G2G2cAWADBH*b[kappa]*Power2(1/(psiexp(1,1) - psiexp(1,2)))*psiexp(-k + nu,1)*psiexp(-Max(0,-j + k),2)*(-(Max(0,-j + k)*(psiexp(1,1) - psiexp(1,2))*psiexp(kappa - nu,2)*psiexp(1 + Max(0,-j + k),1)) + (1 - n + q)*psiexp(1 + kappa - nu,2)*psiexp(1 + Max(0,-j + k),1) + (-1 + kappa - n - nu + q)*psiexp(1 + kappa - nu,1)*psiexp(1 + Max(0,-j + k),2) + (n - q)*psiexp(kappa - nu,2)*psiexp(2 + Max(0,-j + k),1) + (-kappa + n + nu - q)*psiexp(kappa - nu,1)*psiexp(2 + Max(0,-j + k),2)) + vEG2G2G2G2cAWBDBH*b[kappa]*Power2(1/(psiexp(1,1) - psiexp(1,2)))*psiexp(j - k,1)*psiexp(-j + nu - Max(0,-j + k),2)*(Max(0,-j + k)*(psiexp(1,1) - psiexp(1,2))*psiexp(kappa - nu,2)*psiexp(1 + Max(0,-j + k),1) + (-1 + n - q)*psiexp(1 + kappa - nu,2)*psiexp(1 + Max(0,-j + k),1) + (1 - kappa + n + nu - q)*psiexp(1 + kappa - nu,1)*psiexp(1 + Max(0,-j + k),2) + (-n + q)*psiexp(kappa - nu,2)*psiexp(2 + Max(0,-j + k),1) + (kappa - n - nu + q)*psiexp(kappa - nu,1)*psiexp(2 + Max(0,-j + k),2)) - vEG2G2G2G2cAWBDCH*b[kappa]*Power2(1/(psiexp(1,1) - psiexp(1,3)))*psiexp(j - k,1)*psiexp(-j + nu,2)*psiexp(-Max(0,-j + k),3)*(-(Max(0,-j + k)*(psiexp(1,1) - psiexp(1,3))*psiexp(kappa - nu,3)*psiexp(1 + Max(0,-j + k),1)) + (1 - n + q)*psiexp(1 + kappa - nu,3)*psiexp(1 + Max(0,-j + k),1) + (-1 + kappa - n - nu + q)*psiexp(1 + kappa - nu,1)*psiexp(1 + Max(0,-j + k),3) + (n - q)*psiexp(kappa - nu,3)*psiexp(2 + Max(0,-j + k),1) + (-kappa + n + nu - q)*psiexp(kappa - nu,1)*psiexp(2 + Max(0,-j + k),3));
					}
					sumInuP1 = sumIkappa;
				}
				sumInu += sumInuP1*b[nu];
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += MMjkP1/(n - q);
	}
}

for (int j = 1; j <= -2 + q; j++)
{
	for (int k = 1; k <= -3 + q; k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = 2 + j; nu <= Min(q,-1 + j - k + q); nu++)
			{
				scalar sumInuP1;
				{
					scalar sumIkappa = 0;
					for (int kappa = Max(2 + j,1 - j + k + nu); kappa <= q; kappa++)
					{
						sumIkappa += -(vEG2G2G2G2c*(j - kappa + 2*n - 2*q - Max(0,kappa - nu))*(1 + j - kappa + Max(0,kappa - nu))*b[kappa])/2. - (vEG2G2G2G2cAWAD*(j - kappa + 2*n - 2*q - Max(0,kappa - nu))*(1 + j - kappa + Max(0,kappa - nu))*b[kappa]*psiexp(-k + kappa,1))/2. - (vEG2G2G2G2cADAH*(j - kappa + 2*n - 2*q - Max(0,kappa - nu))*(1 + j - kappa + Max(0,kappa - nu))*b[kappa]*psiexp(-j + nu,1))/2. - vEG2G2G2G2cAD*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-Max(0,kappa - nu),1)*((n - q)*psiexp(-j + kappa,1) + Max(0,kappa - nu)*(-1 + psiexp(1,1))*psiexp(-j + kappa,1) + (1 - n + q)*psiexp(1 - j + kappa,1) + (-1 - j + kappa - n + q)*psiexp(1 + Max(0,kappa - nu),1) + (j - kappa + n - q)*psiexp(2 + Max(0,kappa - nu),1)) + vEG2G2G2G2cAW*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(j - k,1)*((-j + kappa - n + q)*psiexp(-j + kappa,1) + (1 + j - kappa + n - q)*psiexp(1 - j + kappa,1) + (-1 + n - q)*psiexp(1 + Max(0,kappa - nu),1) + Max(0,kappa - nu)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,kappa - nu),1) + (-n + q)*psiexp(2 + Max(0,kappa - nu),1)) + vEG2G2G2G2cAWBDBH*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(j - k,1)*psiexp(-j + nu,2)*((-j + kappa - n + q)*psiexp(-j + kappa,1) + (1 + j - kappa + n - q)*psiexp(1 - j + kappa,1) + (-1 + n - q)*psiexp(1 + Max(0,kappa - nu),1) + Max(0,kappa - nu)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,kappa - nu),1) + (-n + q)*psiexp(2 + Max(0,kappa - nu),1)) + vEG2G2G2G2cAWADAH*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-k + nu,1)*((-j + kappa - n + q)*psiexp(-j + kappa,1) + (1 + j - kappa + n - q)*psiexp(1 - j + kappa,1) + (-1 + n - q)*psiexp(1 + Max(0,kappa - nu),1) + Max(0,kappa - nu)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,kappa - nu),1) + (-n + q)*psiexp(2 + Max(0,kappa - nu),1)) + vEG2G2G2G2cAH*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-kappa + nu,1)*((-j + kappa - n + q)*psiexp(-j + kappa,1) + (1 + j - kappa + n - q)*psiexp(1 - j + kappa,1) + (-1 + n - q)*psiexp(1 + Max(0,kappa - nu),1) + Max(0,kappa - nu)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,kappa - nu),1) + (-n + q)*psiexp(2 + Max(0,kappa - nu),1)) + vEG2G2G2G2cAWADBH*b[kappa]*Power2(1/(-1 + psiexp(1,2)))*psiexp(-k + kappa,1)*psiexp(-kappa + nu,2)*((-j + kappa - n + q)*psiexp(-j + kappa,2) + (1 + j - kappa + n - q)*psiexp(1 - j + kappa,2) + (-1 + n - q)*psiexp(1 + Max(0,kappa - nu),2) + Max(0,kappa - nu)*(-1 + psiexp(1,2))*psiexp(1 + Max(0,kappa - nu),2) + (-n + q)*psiexp(2 + Max(0,kappa - nu),2)) - vEG2G2G2G2cADBH*b[kappa]*Power2(1/(psiexp(1,1) - psiexp(1,2)))*psiexp(-kappa + nu,2)*psiexp(-Max(0,kappa - nu),1)*((-1 - j + kappa - n + q)*psiexp(1 - j + kappa,2)*psiexp(1 + Max(0,kappa - nu),1) + Max(0,kappa - nu)*(psiexp(1,1) - psiexp(1,2))*psiexp(-j + kappa,1)*psiexp(1 + Max(0,kappa - nu),2) + (1 - n + q)*psiexp(1 - j + kappa,1)*psiexp(1 + Max(0,kappa - nu),2) + (j - kappa + n - q)*psiexp(-j + kappa,2)*psiexp(2 + Max(0,kappa - nu),1) + (n - q)*psiexp(-j + kappa,1)*psiexp(2 + Max(0,kappa - nu),2)) - vEG2G2G2G2cAWBD*b[kappa]*Power2(1/(psiexp(1,1) - psiexp(1,2)))*psiexp(j - k,1)*psiexp(-Max(0,kappa - nu),2)*(-(Max(0,kappa - nu)*(psiexp(1,1) - psiexp(1,2))*psiexp(-j + kappa,2)*psiexp(1 + Max(0,kappa - nu),1)) + (1 - n + q)*psiexp(1 - j + kappa,2)*psiexp(1 + Max(0,kappa - nu),1) + (-1 - j + kappa - n + q)*psiexp(1 - j + kappa,1)*psiexp(1 + Max(0,kappa - nu),2) + (n - q)*psiexp(-j + kappa,2)*psiexp(2 + Max(0,kappa - nu),1) + (j - kappa + n - q)*psiexp(-j + kappa,1)*psiexp(2 + Max(0,kappa - nu),2)) + vEG2G2G2G2cAWAH*b[kappa]*Power2(1/(-1 + psiexp(2,1)))*psiexp(j - k - kappa + nu,1)*((-j + kappa - n + q)*psiexp(2*(-j + kappa),1) + (1 + j - kappa + n - q)*psiexp(2*(1 - j + kappa),1) + (-1 + n - q)*psiexp(2*(1 + Max(0,kappa - nu)),1) + Max(0,kappa - nu)*(-1 + psiexp(2,1))*psiexp(2*(1 + Max(0,kappa - nu)),1) + (-n + q)*psiexp(2*(2 + Max(0,kappa - nu)),1)) + vEG2G2G2G2cAWBDCH*b[kappa]*Power2(1/(psiexp(1,2) - psiexp(1,1)*psiexp(1,3)))*psiexp(-k,1)*psiexp(-j - kappa + nu,3)*psiexp(-j - Max(0,kappa - nu),2)*((1 + j - kappa + n - q)*psiexp(1 + kappa,1)*psiexp(1 + kappa,3)*psiexp(1 + j + Max(0,kappa - nu),2) + Max(0,kappa - nu)*(-(psiexp(1,2)*psiexp(1 + j,1)) + psiexp(1,3)*psiexp(2 + j,1))*psiexp(kappa,2)*psiexp(Max(0,kappa - nu),1)*psiexp(1 + j + Max(0,kappa - nu),3) + (-1 + n - q)*psiexp(1 + kappa,2)*psiexp(1 + j + Max(0,kappa - nu),1)*psiexp(1 + j + Max(0,kappa - nu),3) + (-j + kappa - n + q)*psiexp(kappa,1)*psiexp(kappa,3)*psiexp(2 + j + Max(0,kappa - nu),2) + (-n + q)*psiexp(kappa,2)*psiexp(2 + j + Max(0,kappa - nu),1)*psiexp(2 + j + Max(0,kappa - nu),3)) - vEG2G2G2G2cAWAHBD*b[kappa]*Power2(1/(-psiexp(1,2) + psiexp(2,1)))*psiexp(j - k - kappa + nu,1)*psiexp(-Max(0,kappa - nu),2)*((-1 - j + kappa - n + q)*psiexp(2 - 2*j + 2*kappa,1)*psiexp(1 + Max(0,kappa - nu),2) - Max(0,kappa - nu)*(-psiexp(1,2) + psiexp(2,1))*psiexp(-j + kappa,2)*psiexp(2*(1 + Max(0,kappa - nu)),1) + (1 - n + q)*psiexp(1 - j + kappa,2)*psiexp(2*(1 + Max(0,kappa - nu)),1) + (j - kappa + n - q)*psiexp(2*(-j + kappa),1)*psiexp(2 + Max(0,kappa - nu),2) + (n - q)*psiexp(-j + kappa,2)*psiexp(4 + 2*Max(0,kappa - nu),1));
					}
					sumInuP1 = sumIkappa;
				}
				sumInu += sumInuP1*b[nu];
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += MMjkP1/(n - q);
	}
}

for (int j = 1; j <= -2 + q; j++)
{
	for (int k = Max(1,2 + 2*j - q); k <= -2 + q; k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = Max(2 + j,2 + 2*j - k); nu <= q; nu++)
			{
				scalar sumInuP1;
				{
					scalar sumIkappa = 0;
					for (int kappa = Max(2 + j,2 + k); kappa <= Min(-j + k + nu,q); kappa++)
					{
						sumIkappa += -(vEG2G2G2G2c*(j - kappa + 2*n - 2*q - Max(0,-j + k))*(1 + j - kappa + Max(0,-j + k))*b[kappa])/2. - (vEG2G2G2G2cAWAD*(j - kappa + 2*n - 2*q - Max(0,-j + k))*(1 + j - kappa + Max(0,-j + k))*b[kappa]*psiexp(-k + kappa,1))/2. - (vEG2G2G2G2cADAH*(j - kappa + 2*n - 2*q - Max(0,-j + k))*(1 + j - kappa + Max(0,-j + k))*b[kappa]*psiexp(-j + nu,1))/2. - vEG2G2G2G2cAD*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-Max(0,-j + k),1)*((n - q)*psiexp(-j + kappa,1) + Max(0,-j + k)*(-1 + psiexp(1,1))*psiexp(-j + kappa,1) + (1 - n + q)*psiexp(1 - j + kappa,1) + (-1 - j + kappa - n + q)*psiexp(1 + Max(0,-j + k),1) + (j - kappa + n - q)*psiexp(2 + Max(0,-j + k),1)) + vEG2G2G2G2cAW*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(j - k,1)*((-j + kappa - n + q)*psiexp(-j + kappa,1) + (1 + j - kappa + n - q)*psiexp(1 - j + kappa,1) + (-1 + n - q)*psiexp(1 + Max(0,-j + k),1) + Max(0,-j + k)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,-j + k),1) + (-n + q)*psiexp(2 + Max(0,-j + k),1)) + vEG2G2G2G2cAWBDBH*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(j - k,1)*psiexp(-j + nu,2)*((-j + kappa - n + q)*psiexp(-j + kappa,1) + (1 + j - kappa + n - q)*psiexp(1 - j + kappa,1) + (-1 + n - q)*psiexp(1 + Max(0,-j + k),1) + Max(0,-j + k)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,-j + k),1) + (-n + q)*psiexp(2 + Max(0,-j + k),1)) + vEG2G2G2G2cAWADAH*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-k + nu,1)*((-j + kappa - n + q)*psiexp(-j + kappa,1) + (1 + j - kappa + n - q)*psiexp(1 - j + kappa,1) + (-1 + n - q)*psiexp(1 + Max(0,-j + k),1) + Max(0,-j + k)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,-j + k),1) + (-n + q)*psiexp(2 + Max(0,-j + k),1)) + vEG2G2G2G2cAH*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-kappa + nu,1)*((-j + kappa - n + q)*psiexp(-j + kappa,1) + (1 + j - kappa + n - q)*psiexp(1 - j + kappa,1) + (-1 + n - q)*psiexp(1 + Max(0,-j + k),1) + Max(0,-j + k)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,-j + k),1) + (-n + q)*psiexp(2 + Max(0,-j + k),1)) + vEG2G2G2G2cAWADBH*b[kappa]*Power2(1/(-1 + psiexp(1,2)))*psiexp(-k + kappa,1)*psiexp(-kappa + nu,2)*((-j + kappa - n + q)*psiexp(-j + kappa,2) + (1 + j - kappa + n - q)*psiexp(1 - j + kappa,2) + (-1 + n - q)*psiexp(1 + Max(0,-j + k),2) + Max(0,-j + k)*(-1 + psiexp(1,2))*psiexp(1 + Max(0,-j + k),2) + (-n + q)*psiexp(2 + Max(0,-j + k),2)) - vEG2G2G2G2cADBH*b[kappa]*Power2(1/(psiexp(1,1) - psiexp(1,2)))*psiexp(-kappa + nu,2)*psiexp(-Max(0,-j + k),1)*((-1 - j + kappa - n + q)*psiexp(1 - j + kappa,2)*psiexp(1 + Max(0,-j + k),1) + Max(0,-j + k)*(psiexp(1,1) - psiexp(1,2))*psiexp(-j + kappa,1)*psiexp(1 + Max(0,-j + k),2) + (1 - n + q)*psiexp(1 - j + kappa,1)*psiexp(1 + Max(0,-j + k),2) + (j - kappa + n - q)*psiexp(-j + kappa,2)*psiexp(2 + Max(0,-j + k),1) + (n - q)*psiexp(-j + kappa,1)*psiexp(2 + Max(0,-j + k),2)) - vEG2G2G2G2cAWBD*b[kappa]*Power2(1/(psiexp(1,1) - psiexp(1,2)))*psiexp(j - k,1)*psiexp(-Max(0,-j + k),2)*(-(Max(0,-j + k)*(psiexp(1,1) - psiexp(1,2))*psiexp(-j + kappa,2)*psiexp(1 + Max(0,-j + k),1)) + (1 - n + q)*psiexp(1 - j + kappa,2)*psiexp(1 + Max(0,-j + k),1) + (-1 - j + kappa - n + q)*psiexp(1 - j + kappa,1)*psiexp(1 + Max(0,-j + k),2) + (n - q)*psiexp(-j + kappa,2)*psiexp(2 + Max(0,-j + k),1) + (j - kappa + n - q)*psiexp(-j + kappa,1)*psiexp(2 + Max(0,-j + k),2)) + vEG2G2G2G2cAWAH*b[kappa]*Power2(1/(-1 + psiexp(2,1)))*psiexp(j - k - kappa + nu,1)*((-j + kappa - n + q)*psiexp(2*(-j + kappa),1) + (1 + j - kappa + n - q)*psiexp(2*(1 - j + kappa),1) + (-1 + n - q)*psiexp(2*(1 + Max(0,-j + k)),1) + Max(0,-j + k)*(-1 + psiexp(2,1))*psiexp(2*(1 + Max(0,-j + k)),1) + (-n + q)*psiexp(2*(2 + Max(0,-j + k)),1)) + vEG2G2G2G2cAWBDCH*b[kappa]*Power2(1/(psiexp(1,2) - psiexp(1,1)*psiexp(1,3)))*psiexp(-k,1)*psiexp(-j - kappa + nu,3)*psiexp(-j - Max(0,-j + k),2)*((1 + j - kappa + n - q)*psiexp(1 + kappa,1)*psiexp(1 + kappa,3)*psiexp(1 + j + Max(0,-j + k),2) + Max(0,-j + k)*(-(psiexp(1,2)*psiexp(1 + j,1)) + psiexp(1,3)*psiexp(2 + j,1))*psiexp(kappa,2)*psiexp(Max(0,-j + k),1)*psiexp(1 + j + Max(0,-j + k),3) + (-1 + n - q)*psiexp(1 + kappa,2)*psiexp(1 + j + Max(0,-j + k),1)*psiexp(1 + j + Max(0,-j + k),3) + (-j + kappa - n + q)*psiexp(kappa,1)*psiexp(kappa,3)*psiexp(2 + j + Max(0,-j + k),2) + (-n + q)*psiexp(kappa,2)*psiexp(2 + j + Max(0,-j + k),1)*psiexp(2 + j + Max(0,-j + k),3)) - vEG2G2G2G2cAWAHBD*b[kappa]*Power2(1/(-psiexp(1,2) + psiexp(2,1)))*psiexp(j - k - kappa + nu,1)*psiexp(-Max(0,-j + k),2)*((-1 - j + kappa - n + q)*psiexp(2 - 2*j + 2*kappa,1)*psiexp(1 + Max(0,-j + k),2) - Max(0,-j + k)*(-psiexp(1,2) + psiexp(2,1))*psiexp(-j + kappa,2)*psiexp(2*(1 + Max(0,-j + k)),1) + (1 - n + q)*psiexp(1 - j + kappa,2)*psiexp(2*(1 + Max(0,-j + k)),1) + (j - kappa + n - q)*psiexp(2*(-j + kappa),1)*psiexp(2 + Max(0,-j + k),2) + (n - q)*psiexp(-j + kappa,2)*psiexp(4 + 2*Max(0,-j + k),1));
					}
					sumInuP1 = sumIkappa;
				}
				sumInu += sumInuP1*b[nu];
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += MMjkP1/(n - q);
	}
}

for (int j = 1; j <= -1 + q; j++)
{
	for (int k = 1; k <= Min(-3 + j + n - q,-1 + q); k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = 1 + j; nu <= q; nu++)
			{
				scalar sumInuP1;
				{
					scalar sumIkappa = 0;
					for (int kappa = 1 + k; kappa <= Min(-2 + j + n - q,q); kappa++)
					{
						sumIkappa += (vEG2G2G2G2c*b[kappa]*(-n + q - 2*n*q + (1 - 2*n + 2*q)*Max(0,-j + kappa) + Power2(n) + Power2(q) + Power2(Max(0,-j + kappa))))/2. + (vEG2G2G2G2cAW*b[kappa]*(-n + q - 2*n*q + (1 - 2*n + 2*q)*Max(0,-j + kappa) + Power2(n) + Power2(q) + Power2(Max(0,-j + kappa)))*psiexp(-k + kappa,1))/2. + (vEG2G2G2G2cAH*b[kappa]*(-n + q - 2*n*q + (1 - 2*n + 2*q)*Max(0,-j + kappa) + Power2(n) + Power2(q) + Power2(Max(0,-j + kappa)))*psiexp(-j + nu,1))/2. + (vEG2G2G2G2cAWAH*b[kappa]*(-n + q - 2*n*q + (1 - 2*n + 2*q)*Max(0,-j + kappa) + Power2(n) + Power2(q) + Power2(Max(0,-j + kappa)))*psiexp(-j - k + kappa + nu,1))/2. + vEG2G2G2G2cAWAD*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(j - k,1)*(psiexp(1 + n - q,1) + (-1 + n - q)*psiexp(1 + Max(0,-j + kappa),1) + Max(0,-j + kappa)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,-j + kappa),1) + (-n + q)*psiexp(2 + Max(0,-j + kappa),1)) + vEG2G2G2G2cAD*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(j - kappa,1)*(psiexp(1 + n - q,1) + (-1 + n - q)*psiexp(1 + Max(0,-j + kappa),1) + Max(0,-j + kappa)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,-j + kappa),1) + (-n + q)*psiexp(2 + Max(0,-j + kappa),1)) + vEG2G2G2G2cAWADBH*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(j - k,1)*psiexp(-j + nu,2)*(psiexp(1 + n - q,1) + (-1 + n - q)*psiexp(1 + Max(0,-j + kappa),1) + Max(0,-j + kappa)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,-j + kappa),1) + (-n + q)*psiexp(2 + Max(0,-j + kappa),1)) + vEG2G2G2G2cADBH*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(j - kappa,1)*psiexp(-j + nu,2)*(psiexp(1 + n - q,1) + (-1 + n - q)*psiexp(1 + Max(0,-j + kappa),1) + Max(0,-j + kappa)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,-j + kappa),1) + (-n + q)*psiexp(2 + Max(0,-j + kappa),1)) + vEG2G2G2G2cAWADAH*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-k + nu,1)*(psiexp(1 + n - q,1) + (-1 + n - q)*psiexp(1 + Max(0,-j + kappa),1) + Max(0,-j + kappa)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,-j + kappa),1) + (-n + q)*psiexp(2 + Max(0,-j + kappa),1)) + vEG2G2G2G2cADAH*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-kappa + nu,1)*(psiexp(1 + n - q,1) + (-1 + n - q)*psiexp(1 + Max(0,-j + kappa),1) + Max(0,-j + kappa)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,-j + kappa),1) + (-n + q)*psiexp(2 + Max(0,-j + kappa),1)) + vEG2G2G2G2cAWBD*b[kappa]*Power2(1/(-1 + psiexp(1,2)))*psiexp(j - kappa,2)*psiexp(-k + kappa,1)*(psiexp(1 + n - q,2) + (-1 + n - q)*psiexp(1 + Max(0,-j + kappa),2) + Max(0,-j + kappa)*(-1 + psiexp(1,2))*psiexp(1 + Max(0,-j + kappa),2) + (-n + q)*psiexp(2 + Max(0,-j + kappa),2)) + vEG2G2G2G2cAWBDCH*b[kappa]*Power2(1/(-1 + psiexp(1,2)))*psiexp(j - kappa,2)*psiexp(-k + kappa,1)*psiexp(-j + nu,3)*(psiexp(1 + n - q,2) + (-1 + n - q)*psiexp(1 + Max(0,-j + kappa),2) + Max(0,-j + kappa)*(-1 + psiexp(1,2))*psiexp(1 + Max(0,-j + kappa),2) + (-n + q)*psiexp(2 + Max(0,-j + kappa),2)) + vEG2G2G2G2cAWBDBH*b[kappa]*Power2(1/(-1 + psiexp(1,2)))*psiexp(-k + kappa,1)*psiexp(-kappa + nu,2)*(psiexp(1 + n - q,2) + (-1 + n - q)*psiexp(1 + Max(0,-j + kappa),2) + Max(0,-j + kappa)*(-1 + psiexp(1,2))*psiexp(1 + Max(0,-j + kappa),2) + (-n + q)*psiexp(2 + Max(0,-j + kappa),2)) + vEG2G2G2G2cAWAHBD*b[kappa]*Power2(1/(-1 + psiexp(1,2)))*psiexp(j - kappa,2)*psiexp(-j - k + kappa + nu,1)*(psiexp(1 + n - q,2) + (-1 + n - q)*psiexp(1 + Max(0,-j + kappa),2) + Max(0,-j + kappa)*(-1 + psiexp(1,2))*psiexp(1 + Max(0,-j + kappa),2) + (-n + q)*psiexp(2 + Max(0,-j + kappa),2));
					}
					sumInuP1 = sumIkappa;
				}
				sumInu += sumInuP1*b[nu];
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += MMjkP1/(n - q);
	}
}

for (int j = 1; j <= -1 + q; j++)
{
	for (int k = 1; k <= Min(-2 + j + n - q,q); k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = 1 + j; nu <= q; nu++)
			{
				scalar sumInuP1;
				{
					scalar sumIkappa = 0;
					for (int kappa = 0; kappa <= -1 + k; kappa++)
					{
						sumIkappa += (vEG2G2G2G2c*b[kappa]*(-n + q - 2*n*q + (1 - 2*n + 2*q)*Max(0,-j + k) + Power2(n) + Power2(q) + Power2(Max(0,-j + k))))/2. + (vEG2G2G2G2cAW*b[kappa]*(-n + q - 2*n*q + (1 - 2*n + 2*q)*Max(0,-j + k) + Power2(n) + Power2(q) + Power2(Max(0,-j + k)))*psiexp(k - kappa,1))/2. + (vEG2G2G2G2cAH*b[kappa]*(-n + q - 2*n*q + (1 - 2*n + 2*q)*Max(0,-j + k) + Power2(n) + Power2(q) + Power2(Max(0,-j + k)))*psiexp(-j + nu,1))/2. + (vEG2G2G2G2cAWAH*b[kappa]*(-n + q - 2*n*q + (1 - 2*n + 2*q)*Max(0,-j + k) + Power2(n) + Power2(q) + Power2(Max(0,-j + k)))*psiexp(-j + k - kappa + nu,1))/2. + vEG2G2G2G2cAD*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(j - k,1)*(psiexp(1 + n - q,1) + (-1 + n - q)*psiexp(1 + Max(0,-j + k),1) + Max(0,-j + k)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,-j + k),1) + (-n + q)*psiexp(2 + Max(0,-j + k),1)) + vEG2G2G2G2cAWAD*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(j - kappa,1)*(psiexp(1 + n - q,1) + (-1 + n - q)*psiexp(1 + Max(0,-j + k),1) + Max(0,-j + k)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,-j + k),1) + (-n + q)*psiexp(2 + Max(0,-j + k),1)) + vEG2G2G2G2cADBH*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(j - k,1)*psiexp(-j + nu,2)*(psiexp(1 + n - q,1) + (-1 + n - q)*psiexp(1 + Max(0,-j + k),1) + Max(0,-j + k)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,-j + k),1) + (-n + q)*psiexp(2 + Max(0,-j + k),1)) + vEG2G2G2G2cAWADBH*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(j - kappa,1)*psiexp(-j + nu,2)*(psiexp(1 + n - q,1) + (-1 + n - q)*psiexp(1 + Max(0,-j + k),1) + Max(0,-j + k)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,-j + k),1) + (-n + q)*psiexp(2 + Max(0,-j + k),1)) + vEG2G2G2G2cADAH*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-k + nu,1)*(psiexp(1 + n - q,1) + (-1 + n - q)*psiexp(1 + Max(0,-j + k),1) + Max(0,-j + k)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,-j + k),1) + (-n + q)*psiexp(2 + Max(0,-j + k),1)) + vEG2G2G2G2cAWADAH*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-kappa + nu,1)*(psiexp(1 + n - q,1) + (-1 + n - q)*psiexp(1 + Max(0,-j + k),1) + Max(0,-j + k)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,-j + k),1) + (-n + q)*psiexp(2 + Max(0,-j + k),1)) + vEG2G2G2G2cAWBD*b[kappa]*Power2(1/(-1 + psiexp(1,2)))*psiexp(j - k,2)*psiexp(k - kappa,1)*(psiexp(1 + n - q,2) + (-1 + n - q)*psiexp(1 + Max(0,-j + k),2) + Max(0,-j + k)*(-1 + psiexp(1,2))*psiexp(1 + Max(0,-j + k),2) + (-n + q)*psiexp(2 + Max(0,-j + k),2)) + vEG2G2G2G2cAWBDCH*b[kappa]*Power2(1/(-1 + psiexp(1,2)))*psiexp(j - k,2)*psiexp(k - kappa,1)*psiexp(-j + nu,3)*(psiexp(1 + n - q,2) + (-1 + n - q)*psiexp(1 + Max(0,-j + k),2) + Max(0,-j + k)*(-1 + psiexp(1,2))*psiexp(1 + Max(0,-j + k),2) + (-n + q)*psiexp(2 + Max(0,-j + k),2)) + vEG2G2G2G2cAWBDBH*b[kappa]*Power2(1/(-1 + psiexp(1,2)))*psiexp(k - kappa,1)*psiexp(-k + nu,2)*(psiexp(1 + n - q,2) + (-1 + n - q)*psiexp(1 + Max(0,-j + k),2) + Max(0,-j + k)*(-1 + psiexp(1,2))*psiexp(1 + Max(0,-j + k),2) + (-n + q)*psiexp(2 + Max(0,-j + k),2)) + vEG2G2G2G2cAWAHBD*b[kappa]*Power2(1/(-1 + psiexp(1,2)))*psiexp(j - k,2)*psiexp(-j + k - kappa + nu,1)*(psiexp(1 + n - q,2) + (-1 + n - q)*psiexp(1 + Max(0,-j + k),2) + Max(0,-j + k)*(-1 + psiexp(1,2))*psiexp(1 + Max(0,-j + k),2) + (-n + q)*psiexp(2 + Max(0,-j + k),2));
					}
					sumInuP1 = sumIkappa;
				}
				sumInu += sumInuP1*b[nu];
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += MMjkP1/(n - q);
	}
}

for (int j = 1; j <= -2 + q; j++)
{
	for (int k = 1; k <= -2 + q; k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = 1 + j; nu <= Min(-1 + q,-1 + j - k + q); nu++)
			{
				scalar sumInuP1;
				{
					scalar sumIkappa = 0;
					for (int kappa = Max(1 + nu,1 - j + k + nu); kappa <= q; kappa++)
					{
						sumIkappa += (-kappa + n + nu - q)*vEG2G2G4c*b[kappa] + (-kappa + n + nu - q)*vEG2G2G4cAWAD*b[kappa]*psiexp(-k + kappa,1) - (kappa - n - nu + q)*vEG2G2G4cAW*b[kappa]*psiexp(j - k + kappa - nu,1) + (-kappa + n + nu - q)*vEG2G2G4cAD*b[kappa]*psiexp(-j + nu,1) - (kappa - n - nu + q)*vEG2G2G4cAWBD*b[kappa]*psiexp(j - k + kappa - nu,1)*psiexp(-j + nu,2);
					}
					sumInuP1 = sumIkappa;
				}
				sumInu += sumInuP1*b[nu];
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += MMjkP1/(n - q);
	}
}

for (int j = 1; j <= -1 + q; j++)
{
	for (int k = 1; k <= -1 + q; k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = 1 + j; nu <= q; nu++)
			{
				scalar sumInuP1;
				{
					scalar sumIkappa = 0;
					for (int kappa = Max(1 + j,1 + k); kappa <= q; kappa++)
					{
						sumIkappa += (j - kappa + n - q)*vEG2G4G2c*b[kappa] + (j - kappa + n - q)*vEG2G4G2cAW*b[kappa]*psiexp(-k + kappa,1) + (j - kappa + n - q)*vEG2G4G2cAH*b[kappa]*psiexp(-j + nu,1) + (j - kappa + n - q)*vEG2G4G2cBH*b[kappa]*psiexp(-j + nu,2) + (j - kappa + n - q)*vEG2G4G2cAWBH*b[kappa]*psiexp(-k + kappa,1)*psiexp(-j + nu,2) + (j - kappa + n - q)*vEG2G4G2cAWCH*b[kappa]*psiexp(-k + kappa,1)*psiexp(-j + nu,3) + (j - kappa + n - q)*vEG2G4G2cAWAH*b[kappa]*psiexp(-j - k + kappa + nu,1);
					}
					sumInuP1 = sumIkappa;
				}
				sumInu += sumInuP1*b[nu];
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += MMjkP1/(n - q);
	}
}

for (int j = 1; j <= -1 + q; j++)
{
	for (int k = 1; k <= Min(-2 + j + n - q,q); k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = 1 + j; nu <= q; nu++)
			{
				sumInu += (vEG4G2G2c*b[nu]*(-n + q - 2*n*q + (1 - 2*n + 2*q)*Max(0,-j + k) + Power2(n) + Power2(q) + Power2(Max(0,-j + k))))/2. + (vEG4G2G2cAH*b[nu]*(-n + q - 2*n*q + (1 - 2*n + 2*q)*Max(0,-j + k) + Power2(n) + Power2(q) + Power2(Max(0,-j + k)))*psiexp(-j + nu,1))/2. + vEG4G2G2cAD*b[nu]*Power2(1/(-1 + psiexp(1,1)))*psiexp(j - k,1)*(psiexp(1 + n - q,1) + (-1 + n - q)*psiexp(1 + Max(0,-j + k),1) + Max(0,-j + k)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,-j + k),1) + (-n + q)*psiexp(2 + Max(0,-j + k),1)) + vEG4G2G2cADBH*b[nu]*Power2(1/(-1 + psiexp(1,1)))*psiexp(j - k,1)*psiexp(-j + nu,2)*(psiexp(1 + n - q,1) + (-1 + n - q)*psiexp(1 + Max(0,-j + k),1) + Max(0,-j + k)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,-j + k),1) + (-n + q)*psiexp(2 + Max(0,-j + k),1)) + vEG4G2G2cADAH*b[nu]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-k + nu,1)*(psiexp(1 + n - q,1) + (-1 + n - q)*psiexp(1 + Max(0,-j + k),1) + Max(0,-j + k)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,-j + k),1) + (-n + q)*psiexp(2 + Max(0,-j + k),1)) + vEG4G2G2cBD*b[nu]*Power2(1/(-1 + psiexp(1,2)))*psiexp(j - k,2)*(psiexp(1 + n - q,2) + (-1 + n - q)*psiexp(1 + Max(0,-j + k),2) + Max(0,-j + k)*(-1 + psiexp(1,2))*psiexp(1 + Max(0,-j + k),2) + (-n + q)*psiexp(2 + Max(0,-j + k),2)) + vEG4G2G2cAHBD*b[nu]*Power2(1/(-1 + psiexp(1,2)))*psiexp(j - k,2)*psiexp(-j + nu,1)*(psiexp(1 + n - q,2) + (-1 + n - q)*psiexp(1 + Max(0,-j + k),2) + Max(0,-j + k)*(-1 + psiexp(1,2))*psiexp(1 + Max(0,-j + k),2) + (-n + q)*psiexp(2 + Max(0,-j + k),2)) + vEG4G2G2cBDCH*b[nu]*Power2(1/(-1 + psiexp(1,2)))*psiexp(j - k,2)*psiexp(-j + nu,3)*(psiexp(1 + n - q,2) + (-1 + n - q)*psiexp(1 + Max(0,-j + k),2) + Max(0,-j + k)*(-1 + psiexp(1,2))*psiexp(1 + Max(0,-j + k),2) + (-n + q)*psiexp(2 + Max(0,-j + k),2)) + vEG4G2G2cBDBH*b[nu]*Power2(1/(-1 + psiexp(1,2)))*psiexp(-k + nu,2)*(psiexp(1 + n - q,2) + (-1 + n - q)*psiexp(1 + Max(0,-j + k),2) + Max(0,-j + k)*(-1 + psiexp(1,2))*psiexp(1 + Max(0,-j + k),2) + (-n + q)*psiexp(2 + Max(0,-j + k),2));
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += (MMjkP1*b[k])/(n - q);
	}
}

for (int j = 1; j <= -3 + q; j++)
{
	for (int k = Max(1,2 + j); k <= -1 + q; k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = 1 + j; nu <= Min(-1 + k,q); nu++)
			{
				scalar sumInuP1;
				{
					scalar sumIkappa = 0;
					for (int kappa = 1 + k; kappa <= q; kappa++)
					{
						sumIkappa += (-k + n + nu - q)*vEG2G4G2c*b[kappa] - (k - n - nu + q)*vEG2G4G2cAH*b[kappa]*psiexp(-k + kappa,1) - (k - n - nu + q)*vEG2G4G2cBH*b[kappa]*psiexp(-k + kappa,2) + (-k + n + nu - q)*vEG2G4G2cAW*b[kappa]*psiexp(-j + nu,1) + (-k + n + nu - q)*vEG2G4G2cAWBH*b[kappa]*psiexp(-k + kappa,2)*psiexp(-j + nu,1) + (-k + n + nu - q)*vEG2G4G2cAWCH*b[kappa]*psiexp(-k + kappa,3)*psiexp(-j + nu,1) - (k - n - nu + q)*vEG2G4G2cAWAH*b[kappa]*psiexp(-j - k + kappa + nu,1);
					}
					sumInuP1 = sumIkappa;
				}
				sumInu += sumInuP1*b[nu];
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += MMjkP1/(n - q);
	}
}

for (int j = 1; j <= -3 + q; j++)
{
	for (int k = Max(1,3 + j); k <= q; k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = 2 + j; nu <= Min(-1 + k,q); nu++)
			{
				scalar sumInuP1;
				{
					scalar sumIkappa = 0;
					for (int kappa = Max(0,1 + j + k - nu); kappa <= -1 + k; kappa++)
					{
						sumIkappa += (-k + n + nu - q)*vEG2G2G4c*b[kappa] - (k - n - nu + q)*vEG2G2G4cAD*b[kappa]*psiexp(k - kappa,1) + (-k + n + nu - q)*vEG2G2G4cAWAD*b[kappa]*psiexp(-j + nu,1) - (k - n - nu + q)*vEG2G2G4cAW*b[kappa]*psiexp(-j - k + kappa + nu,1) - (k - n - nu + q)*vEG2G2G4cAWBD*b[kappa]*psiexp(k - kappa,2)*psiexp(-j - k + kappa + nu,1);
					}
					sumInuP1 = sumIkappa;
				}
				sumInu += sumInuP1*b[nu];
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += MMjkP1/(n - q);
	}
}

for (int j = 1; j <= -2 + q; j++)
{
	for (int k = Max(2,2 + j); k <= q; k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = 1 + j; nu <= Min3(-1 + k,-1 + j + k,q); nu++)
			{
				scalar sumInuP1;
				{
					scalar sumIkappa = 0;
					for (int kappa = 0; kappa <= -1 + j + k - nu; kappa++)
					{
						sumIkappa += (-k + n + nu - q)*vEG2G2G4c*b[kappa] - (k - n - nu + q)*vEG2G2G4cAWAD*b[kappa]*psiexp(k - kappa,1) - (k - n - nu + q)*vEG2G2G4cAW*b[kappa]*psiexp(j + k - kappa - nu,1) + (-k + n + nu - q)*vEG2G2G4cAD*b[kappa]*psiexp(-j + nu,1) - (k - n - nu + q)*vEG2G2G4cAWBD*b[kappa]*psiexp(j + k - kappa - nu,1)*psiexp(-j + nu,2);
					}
					sumInuP1 = sumIkappa;
				}
				sumInu += sumInuP1*b[nu];
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += MMjkP1/(n - q);
	}
}

for (int j = 1; j <= -2 + q; j++)
{
	for (int k = Max(1,2 + j); k <= q; k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = 1 + j; nu <= Min(-1 + k,q); nu++)
			{
				sumInu += (-k + n + nu - q)*vEG2G6c*b[nu] + (-k + n + nu - q)*vEG2G6cAW*b[nu]*psiexp(-j + nu,1);
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += (MMjkP1*b[k])/(n - q);
	}
}

for (int j = 1; j <= -2 + q; j++)
{
	for (int k = Max(1,2 + j); k <= q; k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = 1 + j; nu <= Min3(-1 + k,j + k,q); nu++)
			{
				scalar sumInuP1;
				{
					scalar sumIkappa = 0;
					for (int kappa = Max(0,j + k - nu); kappa <= j + k - nu; kappa++)
					{
						sumIkappa += b[kappa];
					}
					sumInuP1 = sumIkappa;
				}
				sumInu += (-k + n + nu - q)*sumInuP1*b[nu]*(vEG4G4c + vEG4G4cAD*psiexp(-j + nu,1) + vEG4G4cBD*psiexp(-j + nu,2));
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += MMjkP1/(n - q);
	}
}

for (int j = 1; j <= -3 + q; j++)
{
	for (int k = Max(1,1 + j); k <= -2 + q; k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = 1 + j; nu <= Min(q,-1 + j - k + q); nu++)
			{
				scalar sumInuP1;
				{
					scalar sumIkappa = 0;
					for (int kappa = 1 - j + k + nu; kappa <= q; kappa++)
					{
						sumIkappa += (j - k + n - q)*vEG4G2G2c*b[kappa] + (j - k + n - q)*vEG4G2G2cADAH*b[kappa]*psiexp(-k + kappa,1) + (j - k + n - q)*vEG4G2G2cBDBH*b[kappa]*psiexp(-k + kappa,2) + (j - k + n - q)*vEG4G2G2cAH*b[kappa]*psiexp(j - k + kappa - nu,1) + (j - k + n - q)*vEG4G2G2cAD*b[kappa]*psiexp(-j + nu,1) + (j - k + n - q)*vEG4G2G2cADBH*b[kappa]*psiexp(j - k + kappa - nu,2)*psiexp(-j + nu,1) + (j - k + n - q)*vEG4G2G2cBD*b[kappa]*psiexp(-j + nu,2) + (j - k + n - q)*vEG4G2G2cAHBD*b[kappa]*psiexp(j - k + kappa - nu,1)*psiexp(-j + nu,2) + (j - k + n - q)*vEG4G2G2cBDCH*b[kappa]*psiexp(j - k + kappa - nu,3)*psiexp(-j + nu,2);
					}
					sumInuP1 = sumIkappa;
				}
				sumInu += sumInuP1*b[nu];
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += MMjkP1/(n - q);
	}
}

for (int j = 1; j <= -2 + q; j++)
{
	for (int k = Max(1,1 + j); k <= -1 + q; k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = 2 + j; nu <= q; nu++)
			{
				scalar sumInuP1;
				{
					scalar sumIkappa = 0;
					for (int kappa = 1 + k; kappa <= Min(-1 - j + k + nu,q); kappa++)
					{
						sumIkappa += (j - k + n - q)*vEG4G2G2c*b[kappa] + (j - k + n - q)*vEG4G2G2cAD*b[kappa]*psiexp(-k + kappa,1) + (j - k + n - q)*vEG4G2G2cBD*b[kappa]*psiexp(-k + kappa,2) + (j - k + n - q)*vEG4G2G2cADAH*b[kappa]*psiexp(-j + nu,1) + (j - k + n - q)*vEG4G2G2cBDBH*b[kappa]*psiexp(-j + nu,2) + (j - k + n - q)*vEG4G2G2cAH*b[kappa]*psiexp(-j + k - kappa + nu,1) + (j - k + n - q)*vEG4G2G2cAHBD*b[kappa]*psiexp(-k + kappa,2)*psiexp(-j + k - kappa + nu,1) + (j - k + n - q)*vEG4G2G2cADBH*b[kappa]*psiexp(-k + kappa,1)*psiexp(-j + k - kappa + nu,2) + (j - k + n - q)*vEG4G2G2cBDCH*b[kappa]*psiexp(-k + kappa,2)*psiexp(-j + k - kappa + nu,3);
					}
					sumInuP1 = sumIkappa;
				}
				sumInu += sumInuP1*b[nu];
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += MMjkP1/(n - q);
	}
}

for (int j = 1; j <= -1 + q; j++)
{
	for (int k = Max(1,1 + j); k <= q; k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = 1 + j; nu <= q; nu++)
			{
				scalar sumInuP1;
				{
					scalar sumIkappa = 0;
					for (int kappa = 0; kappa <= -1 + k; kappa++)
					{
						sumIkappa += (j - k + n - q)*vEG2G4G2c*b[kappa] + (j - k + n - q)*vEG2G4G2cAW*b[kappa]*psiexp(k - kappa,1) + (j - k + n - q)*vEG2G4G2cAH*b[kappa]*psiexp(-j + nu,1) + (j - k + n - q)*vEG2G4G2cBH*b[kappa]*psiexp(-j + nu,2) + (j - k + n - q)*vEG2G4G2cAWBH*b[kappa]*psiexp(k - kappa,1)*psiexp(-j + nu,2) + (j - k + n - q)*vEG2G4G2cAWCH*b[kappa]*psiexp(k - kappa,1)*psiexp(-j + nu,3) + (j - k + n - q)*vEG2G4G2cAWAH*b[kappa]*psiexp(-j + k - kappa + nu,1);
					}
					sumInuP1 = sumIkappa;
				}
				sumInu += sumInuP1*b[nu];
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += MMjkP1/(n - q);
	}
}

for (int j = 1; j <= -2 + q; j++)
{
	for (int k = Max(1,1 + j); k <= -1 + q; k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = 1 + j; nu <= Min(q,j - k + q); nu++)
			{
				scalar sumInuP1;
				{
					scalar sumIkappa = 0;
					for (int kappa = -j + k + nu; kappa <= Min(-j + k + nu,q); kappa++)
					{
						sumIkappa += b[kappa];
					}
					sumInuP1 = sumIkappa;
				}
				sumInu += sumInuP1*b[nu]*(vEG4G4c + vEG4G4cAD*psiexp(-j + nu,1) + vEG4G4cBD*psiexp(-j + nu,2));
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += (MMjkP1*(j - k + n - q))/(n - q);
	}
}

for (int j = 1; j <= -1 + q; j++)
{
	for (int k = Max(1,1 + j); k <= q; k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = 1 + j; nu <= q; nu++)
			{
				sumInu += vEG6G2c*b[nu] + vEG6G2cAH*b[nu]*psiexp(-j + nu,1) + vEG6G2cBH*b[nu]*psiexp(-j + nu,2) + vEG6G2cCH*b[nu]*psiexp(-j + nu,3);
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += (MMjkP1*(j - k + n - q)*b[k])/(n - q);
	}
}

for (int j = 1; j <= -3 + q; j++)
{
	for (int k = Max(1,2 + j); k <= -1 + q; k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = 0; nu <= -1 + j; nu++)
			{
				scalar sumInuP1;
				{
					scalar sumIkappa = 0;
					for (int kappa = 1 + k; kappa <= q; kappa++)
					{
						sumIkappa += -((1 + j - k)*(j - k + 2*n - 2*q)*vEG2G2G2G2c*b[kappa])/2. + vEG2G2G2G2cAD*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*(-n + q + (-1 + n - q)*psiexp(1,1) + (1 + j - k + n - q)*psiexp(1 + j - k,1) + (-j + k - n + q)*psiexp(2 + j - k,1))*psiexp(-j + k,1) + vEG2G2G2G2cADAH*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*(-n + q + (-1 + n - q)*psiexp(1,1) + (1 + j - k + n - q)*psiexp(1 + j - k,1) + (-j + k - n + q)*psiexp(2 + j - k,1))*psiexp(-j + kappa,1) - ((1 + j - k)*(j - k + 2*n - 2*q)*vEG2G2G2G2cAH*b[kappa]*psiexp(-k + kappa,1))/2. + vEG2G2G2G2cADBH*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*(-n + q + (-1 + n - q)*psiexp(1,1) + (1 + j - k + n - q)*psiexp(1 + j - k,1) + (-j + k - n + q)*psiexp(2 + j - k,1))*psiexp(-j + k,1)*psiexp(-k + kappa,2) - ((1 + j - k)*(j - k + 2*n - 2*q)*vEG2G2G2G2cAW*b[kappa]*psiexp(j - nu,1))/2. + vEG2G2G2G2cAWBD*b[kappa]*Power2(1/(-1 + psiexp(1,2)))*(-n + q + (-1 + n - q)*psiexp(1,2) + (1 + j - k + n - q)*psiexp(1 + j - k,2) + (-j + k - n + q)*psiexp(2 + j - k,2))*psiexp(-j + k,2)*psiexp(j - nu,1) + vEG2G2G2G2cAWBDBH*b[kappa]*Power2(1/(-1 + psiexp(1,2)))*(-n + q + (-1 + n - q)*psiexp(1,2) + (1 + j - k + n - q)*psiexp(1 + j - k,2) + (-j + k - n + q)*psiexp(2 + j - k,2))*psiexp(-j + kappa,2)*psiexp(j - nu,1) + vEG2G2G2G2cAWBDCH*b[kappa]*Power2(1/(-1 + psiexp(1,2)))*(-n + q + (-1 + n - q)*psiexp(1,2) + (1 + j - k + n - q)*psiexp(1 + j - k,2) + (-j + k - n + q)*psiexp(2 + j - k,2))*psiexp(-j + k,2)*psiexp(-k + kappa,3)*psiexp(j - nu,1) + vEG2G2G2G2cAWAD*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*(-n + q + (-1 + n - q)*psiexp(1,1) + (1 + j - k + n - q)*psiexp(1 + j - k,1) + (-j + k - n + q)*psiexp(2 + j - k,1))*psiexp(k - nu,1) + vEG2G2G2G2cAWADBH*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*(-n + q + (-1 + n - q)*psiexp(1,1) + (1 + j - k + n - q)*psiexp(1 + j - k,1) + (-j + k - n + q)*psiexp(2 + j - k,1))*psiexp(-k + kappa,2)*psiexp(k - nu,1) + vEG2G2G2G2cAWADAH*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*(-n + q + (-1 + n - q)*psiexp(1,1) + (1 + j - k + n - q)*psiexp(1 + j - k,1) + (-j + k - n + q)*psiexp(2 + j - k,1))*psiexp(kappa - nu,1) - ((1 + j - k)*(j - k + 2*n - 2*q)*vEG2G2G2G2cAWAH*b[kappa]*psiexp(j - k + kappa - nu,1))/2. + vEG2G2G2G2cAWAHBD*b[kappa]*Power2(1/(-1 + psiexp(1,2)))*(-n + q + (-1 + n - q)*psiexp(1,2) + (1 + j - k + n - q)*psiexp(1 + j - k,2) + (-j + k - n + q)*psiexp(2 + j - k,2))*psiexp(-j + k,2)*psiexp(j - k + kappa - nu,1);
					}
					sumInuP1 = sumIkappa;
				}
				sumInu += sumInuP1*b[nu];
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += MMjkP1/(n - q);
	}
}

for (int j = 1; j <= -3 + q; j++)
{
	for (int k = Max(1,3 + j); k <= q; k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = 0; nu <= -1 + j; nu++)
			{
				scalar sumInuP1;
				{
					scalar sumIkappa = 0;
					for (int kappa = Max(0,2 + j); kappa <= -1 + k; kappa++)
					{
						sumIkappa += -((1 + j - kappa)*(j - kappa + 2*n - 2*q)*vEG2G2G2G2c*b[kappa])/2. + vEG2G2G2G2cADAH*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-j + k,1)*(-n + q + (-1 + n - q)*psiexp(1,1) + (1 + j - kappa + n - q)*psiexp(1 + j - kappa,1) + (-j + kappa - n + q)*psiexp(2 + j - kappa,1)) - ((-1 - j + kappa)*(-j + kappa - 2*n + 2*q)*vEG2G2G2G2cAH*b[kappa]*psiexp(k - kappa,1))/2. + vEG2G2G2G2cAD*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*(-n + q + (-1 + n - q)*psiexp(1,1) + (1 + j - kappa + n - q)*psiexp(1 + j - kappa,1) + (-j + kappa - n + q)*psiexp(2 + j - kappa,1))*psiexp(-j + kappa,1) + vEG2G2G2G2cADBH*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*(-n + q + (-1 + n - q)*psiexp(1,1) + (1 + j - kappa + n - q)*psiexp(1 + j - kappa,1) + (-j + kappa - n + q)*psiexp(2 + j - kappa,1))*psiexp(k - kappa,2)*psiexp(-j + kappa,1) - ((1 + j - kappa)*(j - kappa + 2*n - 2*q)*vEG2G2G2G2cAW*b[kappa]*psiexp(j - nu,1))/2. + vEG2G2G2G2cAWBDBH*b[kappa]*Power2(1/(-1 + psiexp(1,2)))*psiexp(-j + k,2)*(-n + q + (-1 + n - q)*psiexp(1,2) + (1 + j - kappa + n - q)*psiexp(1 + j - kappa,2) + (-j + kappa - n + q)*psiexp(2 + j - kappa,2))*psiexp(j - nu,1) + vEG2G2G2G2cAWBD*b[kappa]*Power2(1/(-1 + psiexp(1,2)))*(-n + q + (-1 + n - q)*psiexp(1,2) + (1 + j - kappa + n - q)*psiexp(1 + j - kappa,2) + (-j + kappa - n + q)*psiexp(2 + j - kappa,2))*psiexp(-j + kappa,2)*psiexp(j - nu,1) + vEG2G2G2G2cAWBDCH*b[kappa]*Power2(1/(-1 + psiexp(1,2)))*(-n + q + (-1 + n - q)*psiexp(1,2) + (1 + j - kappa + n - q)*psiexp(1 + j - kappa,2) + (-j + kappa - n + q)*psiexp(2 + j - kappa,2))*psiexp(k - kappa,3)*psiexp(-j + kappa,2)*psiexp(j - nu,1) + vEG2G2G2G2cAWADAH*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*(-n + q + (-1 + n - q)*psiexp(1,1) + (1 + j - kappa + n - q)*psiexp(1 + j - kappa,1) + (-j + kappa - n + q)*psiexp(2 + j - kappa,1))*psiexp(k - nu,1) - ((1 + j - kappa)*(j - kappa + 2*n - 2*q)*vEG2G2G2G2cAWAH*b[kappa]*psiexp(j + k - kappa - nu,1))/2. + vEG2G2G2G2cAWAHBD*b[kappa]*Power2(1/(-1 + psiexp(1,2)))*(-n + q + (-1 + n - q)*psiexp(1,2) + (1 + j - kappa + n - q)*psiexp(1 + j - kappa,2) + (-j + kappa - n + q)*psiexp(2 + j - kappa,2))*psiexp(-j + kappa,2)*psiexp(j + k - kappa - nu,1) + vEG2G2G2G2cAWAD*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*(-n + q + (-1 + n - q)*psiexp(1,1) + (1 + j - kappa + n - q)*psiexp(1 + j - kappa,1) + (-j + kappa - n + q)*psiexp(2 + j - kappa,1))*psiexp(kappa - nu,1) + vEG2G2G2G2cAWADBH*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*(-n + q + (-1 + n - q)*psiexp(1,1) + (1 + j - kappa + n - q)*psiexp(1 + j - kappa,1) + (-j + kappa - n + q)*psiexp(2 + j - kappa,1))*psiexp(k - kappa,2)*psiexp(kappa - nu,1);
					}
					sumInuP1 = sumIkappa;
				}
				sumInu += sumInuP1*b[nu];
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += MMjkP1/(n - q);
	}
}

for (int j = 2; j <= -2 + q; j++)
{
	for (int k = Max(2,2 + j); k <= q; k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = Max(0,j - k); nu <= -2 + j; nu++)
			{
				scalar sumInuP1;
				{
					scalar sumIkappa = 0;
					for (int kappa = Max(0,2 + nu); kappa <= -j + k + nu; kappa++)
					{
						sumIkappa += -(vEG2G2G2G2c*(-1 + kappa - nu - Max(0,-j + kappa))*(kappa - 2*n - nu + 2*q + Max(0,-j + kappa))*b[kappa])/2. - (vEG2G2G2G2cADAH*(-1 + kappa - nu - Max(0,-j + kappa))*(kappa - 2*n - nu + 2*q + Max(0,-j + kappa))*b[kappa]*psiexp(k - kappa,1))/2. - (vEG2G2G2G2cAWAD*(-kappa + 2*n + nu - 2*q - Max(0,-j + kappa))*(1 - kappa + nu + Max(0,-j + kappa))*b[kappa]*psiexp(j - nu,1))/2. + vEG2G2G2G2cAD*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(j - kappa,1)*((kappa - n - nu + q)*psiexp(kappa - nu,1) + (1 - kappa + n + nu - q)*psiexp(1 + kappa - nu,1) + (-1 + n - q)*psiexp(1 + Max(0,-j + kappa),1) + Max(0,-j + kappa)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,-j + kappa),1) + (-n + q)*psiexp(2 + Max(0,-j + kappa),1)) + vEG2G2G2G2cAWADAH*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(k - kappa - Max(0,-j + kappa),1)*((-n + q)*psiexp(kappa - nu,1) - Max(0,-j + kappa)*(-1 + psiexp(1,1))*psiexp(kappa - nu,1) + (-1 + n - q)*psiexp(1 + kappa - nu,1) + (1 - kappa + n + nu - q)*psiexp(1 + Max(0,-j + kappa),1) + (kappa - n - nu + q)*psiexp(2 + Max(0,-j + kappa),1)) + vEG2G2G2G2cAH*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-j + k - kappa + nu - Max(0,-j + kappa),1)*((-n + q)*psiexp(kappa - nu,1) - Max(0,-j + kappa)*(-1 + psiexp(1,1))*psiexp(kappa - nu,1) + (-1 + n - q)*psiexp(1 + kappa - nu,1) + (1 - kappa + n + nu - q)*psiexp(1 + Max(0,-j + kappa),1) + (kappa - n - nu + q)*psiexp(2 + Max(0,-j + kappa),1)) + vEG2G2G2G2cAW*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-Max(0,-j + kappa),1)*((-n + q)*psiexp(kappa - nu,1) - Max(0,-j + kappa)*(-1 + psiexp(1,1))*psiexp(kappa - nu,1) + (-1 + n - q)*psiexp(1 + kappa - nu,1) + (1 - kappa + n + nu - q)*psiexp(1 + Max(0,-j + kappa),1) + (kappa - n - nu + q)*psiexp(2 + Max(0,-j + kappa),1)) + vEG2G2G2G2cAWBDBH*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(k - kappa,2)*psiexp(-Max(0,-j + kappa),1)*((-n + q)*psiexp(kappa - nu,1) - Max(0,-j + kappa)*(-1 + psiexp(1,1))*psiexp(kappa - nu,1) + (-1 + n - q)*psiexp(1 + kappa - nu,1) + (1 - kappa + n + nu - q)*psiexp(1 + Max(0,-j + kappa),1) + (kappa - n - nu + q)*psiexp(2 + Max(0,-j + kappa),1)) - vEG2G2G2G2cAWADBH*b[kappa]*Power2(1/(-1 + psiexp(1,2)))*psiexp(j - nu,1)*psiexp(-j + k - kappa + nu - Max(0,-j + kappa),2)*((n - q)*psiexp(kappa - nu,2) + Max(0,-j + kappa)*(-1 + psiexp(1,2))*psiexp(kappa - nu,2) + (1 - n + q)*psiexp(1 + kappa - nu,2) + (-1 + kappa - n - nu + q)*psiexp(1 + Max(0,-j + kappa),2) + (-kappa + n + nu - q)*psiexp(2 + Max(0,-j + kappa),2)) + vEG2G2G2G2cAWBD*b[kappa]*Power2(1/(psiexp(1,1) - psiexp(1,2)))*psiexp(j - kappa,2)*psiexp(-Max(0,-j + kappa),1)*((1 - kappa + n + nu - q)*psiexp(1 + kappa - nu,2)*psiexp(1 + Max(0,-j + kappa),1) + Max(0,-j + kappa)*(-psiexp(1,1) + psiexp(1,2))*psiexp(kappa - nu,1)*psiexp(1 + Max(0,-j + kappa),2) + (-1 + n - q)*psiexp(1 + kappa - nu,1)*psiexp(1 + Max(0,-j + kappa),2) + (kappa - n - nu + q)*psiexp(kappa - nu,2)*psiexp(2 + Max(0,-j + kappa),1) + (-n + q)*psiexp(kappa - nu,1)*psiexp(2 + Max(0,-j + kappa),2)) + vEG2G2G2G2cADBH*b[kappa]*Power2(1/(psiexp(1,1) - psiexp(1,2)))*psiexp(j - kappa,1)*psiexp(-j + k - kappa + nu - Max(0,-j + kappa),2)*(Max(0,-j + kappa)*(psiexp(1,1) - psiexp(1,2))*psiexp(kappa - nu,2)*psiexp(1 + Max(0,-j + kappa),1) + (-1 + n - q)*psiexp(1 + kappa - nu,2)*psiexp(1 + Max(0,-j + kappa),1) + (1 - kappa + n + nu - q)*psiexp(1 + kappa - nu,1)*psiexp(1 + Max(0,-j + kappa),2) + (-n + q)*psiexp(kappa - nu,2)*psiexp(2 + Max(0,-j + kappa),1) + (kappa - n - nu + q)*psiexp(kappa - nu,1)*psiexp(2 + Max(0,-j + kappa),2)) + vEG2G2G2G2cAWAH*b[kappa]*Power2(1/(-1 + psiexp(2,1)))*psiexp(-j + k - kappa + nu - 2*Max(0,-j + kappa),1)*((-n + q)*psiexp(2*(kappa - nu),1) - Max(0,-j + kappa)*(-1 + psiexp(2,1))*psiexp(2*(kappa - nu),1) + (-1 + n - q)*psiexp(2*(1 + kappa - nu),1) + (1 - kappa + n + nu - q)*psiexp(2*(1 + Max(0,-j + kappa)),1) + (kappa - n - nu + q)*psiexp(2*(2 + Max(0,-j + kappa)),1)) + vEG2G2G2G2cAWBDCH*b[kappa]*Power2(1/(psiexp(1,2) - psiexp(1,1)*psiexp(1,3)))*psiexp(j - kappa - nu,2)*psiexp(-j + k - kappa - Max(0,-j + kappa),3)*psiexp(-nu - Max(0,-j + kappa),1)*(Max(0,-j + kappa)*(psiexp(1,2) - psiexp(1,1)*psiexp(1,3))*psiexp(kappa,1)*psiexp(kappa,3)*psiexp(1 + nu + Max(0,-j + kappa),2) + (-1 + n - q)*psiexp(1 + kappa,1)*psiexp(1 + kappa,3)*psiexp(1 + nu + Max(0,-j + kappa),2) + (1 - kappa + n + nu - q)*psiexp(1 + kappa,2)*psiexp(1 + nu + Max(0,-j + kappa),1)*psiexp(1 + nu + Max(0,-j + kappa),3) + (-n + q)*psiexp(kappa,1)*psiexp(kappa,3)*psiexp(2 + nu + Max(0,-j + kappa),2) + (kappa - n - nu + q)*psiexp(kappa,2)*psiexp(2 + nu + Max(0,-j + kappa),1)*psiexp(2 + nu + Max(0,-j + kappa),3)) + vEG2G2G2G2cAWAHBD*b[kappa]*Power2(1/(-psiexp(1,2) + psiexp(2,1)))*psiexp(j - kappa,2)*psiexp(-j + k - kappa + nu - 2*Max(0,-j + kappa),1)*(Max(0,-j + kappa)*(psiexp(1,2) - psiexp(2,1))*psiexp(2*kappa - 2*nu,1)*psiexp(1 + Max(0,-j + kappa),2) + (-1 + n - q)*psiexp(2*(1 + kappa - nu),1)*psiexp(1 + Max(0,-j + kappa),2) + (1 - kappa + n + nu - q)*psiexp(1 + kappa - nu,2)*psiexp(2*(1 + Max(0,-j + kappa)),1) + (-n + q)*psiexp(2*(kappa - nu),1)*psiexp(2 + Max(0,-j + kappa),2) + (kappa - n - nu + q)*psiexp(kappa - nu,2)*psiexp(4 + 2*Max(0,-j + kappa),1));
					}
					sumInuP1 = sumIkappa;
				}
				sumInu += sumInuP1*b[nu];
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += MMjkP1/(n - q);
	}
}

for (int j = 3; j <= -2 + q; j++)
{
	for (int k = Max(2,2 + j); k <= q; k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = 0; nu <= -3 + j; nu++)
			{
				scalar sumInuP1;
				{
					scalar sumIkappa = 0;
					for (int kappa = Max(0,1 - j + k + nu); kappa <= -2 + k; kappa++)
					{
						sumIkappa += -(vEG2G2G2G2c*(j - k + 2*n - 2*q - Max(0,-j + kappa))*(1 + j - k + Max(0,-j + kappa))*b[kappa])/2. - (vEG2G2G2G2cADAH*(-1 - j + k - Max(0,-j + kappa))*(-j + k - 2*n + 2*q + Max(0,-j + kappa))*b[kappa]*psiexp(k - kappa,1))/2. - (vEG2G2G2G2cAWAD*(j - k + 2*n - 2*q - Max(0,-j + kappa))*(1 + j - k + Max(0,-j + kappa))*b[kappa]*psiexp(j - nu,1))/2. - vEG2G2G2G2cAW*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(j - k + kappa - nu - Max(0,-j + kappa),1)*((n - q)*psiexp(-j + k,1) + Max(0,-j + kappa)*(-1 + psiexp(1,1))*psiexp(-j + k,1) + (1 - n + q)*psiexp(1 - j + k,1) + (-1 - j + k - n + q)*psiexp(1 + Max(0,-j + kappa),1) + (j - k + n - q)*psiexp(2 + Max(0,-j + kappa),1)) - vEG2G2G2G2cAH*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-Max(0,-j + kappa),1)*((n - q)*psiexp(-j + k,1) + Max(0,-j + kappa)*(-1 + psiexp(1,1))*psiexp(-j + k,1) + (1 - n + q)*psiexp(1 - j + k,1) + (-1 - j + k - n + q)*psiexp(1 + Max(0,-j + kappa),1) + (j - k + n - q)*psiexp(2 + Max(0,-j + kappa),1)) + vEG2G2G2G2cAD*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(j - kappa,1)*((-j + k - n + q)*psiexp(-j + k,1) + (1 + j - k + n - q)*psiexp(1 - j + k,1) + (-1 + n - q)*psiexp(1 + Max(0,-j + kappa),1) + Max(0,-j + kappa)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,-j + kappa),1) + (-n + q)*psiexp(2 + Max(0,-j + kappa),1)) + vEG2G2G2G2cAWADAH*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(j - nu - Max(0,-j + kappa),1)*((-n + q)*psiexp(-j + k,1) - Max(0,-j + kappa)*(-1 + psiexp(1,1))*psiexp(-j + k,1) + (-1 + n - q)*psiexp(1 - j + k,1) + (1 + j - k + n - q)*psiexp(1 + Max(0,-j + kappa),1) + (-j + k - n + q)*psiexp(2 + Max(0,-j + kappa),1)) + vEG2G2G2G2cAWBDBH*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(k - kappa,2)*psiexp(j - k + kappa - nu - Max(0,-j + kappa),1)*((-n + q)*psiexp(-j + k,1) - Max(0,-j + kappa)*(-1 + psiexp(1,1))*psiexp(-j + k,1) + (-1 + n - q)*psiexp(1 - j + k,1) + (1 + j - k + n - q)*psiexp(1 + Max(0,-j + kappa),1) + (-j + k - n + q)*psiexp(2 + Max(0,-j + kappa),1)) + vEG2G2G2G2cAWADBH*b[kappa]*Power2(1/(-1 + psiexp(1,2)))*psiexp(j - nu,1)*psiexp(-Max(0,-j + kappa),2)*((-n + q)*psiexp(-j + k,2) - Max(0,-j + kappa)*(-1 + psiexp(1,2))*psiexp(-j + k,2) + (-1 + n - q)*psiexp(1 - j + k,2) + (1 + j - k + n - q)*psiexp(1 + Max(0,-j + kappa),2) + (-j + k - n + q)*psiexp(2 + Max(0,-j + kappa),2)) - vEG2G2G2G2cADBH*b[kappa]*Power2(1/(psiexp(1,1) - psiexp(1,2)))*psiexp(j - kappa,1)*psiexp(-Max(0,-j + kappa),2)*(-(Max(0,-j + kappa)*(psiexp(1,1) - psiexp(1,2))*psiexp(-j + k,2)*psiexp(1 + Max(0,-j + kappa),1)) + (1 - n + q)*psiexp(1 - j + k,2)*psiexp(1 + Max(0,-j + kappa),1) + (-1 - j + k - n + q)*psiexp(1 - j + k,1)*psiexp(1 + Max(0,-j + kappa),2) + (n - q)*psiexp(-j + k,2)*psiexp(2 + Max(0,-j + kappa),1) + (j - k + n - q)*psiexp(-j + k,1)*psiexp(2 + Max(0,-j + kappa),2)) + vEG2G2G2G2cAWBD*b[kappa]*Power2(1/(psiexp(1,1) - psiexp(1,2)))*psiexp(j - kappa,2)*psiexp(j - k + kappa - nu - Max(0,-j + kappa),1)*((1 + j - k + n - q)*psiexp(1 - j + k,2)*psiexp(1 + Max(0,-j + kappa),1) + Max(0,-j + kappa)*(-psiexp(1,1) + psiexp(1,2))*psiexp(-j + k,1)*psiexp(1 + Max(0,-j + kappa),2) + (-1 + n - q)*psiexp(1 - j + k,1)*psiexp(1 + Max(0,-j + kappa),2) + (-j + k - n + q)*psiexp(-j + k,2)*psiexp(2 + Max(0,-j + kappa),1) + (-n + q)*psiexp(-j + k,1)*psiexp(2 + Max(0,-j + kappa),2)) + vEG2G2G2G2cAWAH*b[kappa]*Power2(1/(-1 + psiexp(2,1)))*psiexp(j - k + kappa - nu - 2*Max(0,-j + kappa),1)*((-n + q)*psiexp(2*(-j + k),1) - Max(0,-j + kappa)*(-1 + psiexp(2,1))*psiexp(2*(-j + k),1) + (-1 + n - q)*psiexp(2*(1 - j + k),1) + (1 + j - k + n - q)*psiexp(2*(1 + Max(0,-j + kappa)),1) + (-j + k - n + q)*psiexp(2*(2 + Max(0,-j + kappa)),1)) - vEG2G2G2G2cAWBDCH*b[kappa]*Power2(1/(psiexp(1,2) - psiexp(1,1)*psiexp(1,3)))*psiexp(-kappa,2)*psiexp(-j - Max(0,-j + kappa),3)*psiexp(-k + kappa - nu - Max(0,-j + kappa),1)*(Max(0,-j + kappa)*(-psiexp(1,2) + psiexp(1,1)*psiexp(1,3))*psiexp(k,1)*psiexp(k,3)*psiexp(1 + j + Max(0,-j + kappa),2) + (1 - n + q)*psiexp(1 + k,1)*psiexp(1 + k,3)*psiexp(1 + j + Max(0,-j + kappa),2) + (-1 - j + k - n + q)*psiexp(1 + k,2)*psiexp(1 + j + Max(0,-j + kappa),1)*psiexp(1 + j + Max(0,-j + kappa),3) + (n - q)*psiexp(k,1)*psiexp(k,3)*psiexp(2 + j + Max(0,-j + kappa),2) + (j - k + n - q)*psiexp(k,2)*psiexp(2 + j + Max(0,-j + kappa),1)*psiexp(2 + j + Max(0,-j + kappa),3)) - vEG2G2G2G2cAWAHBD*b[kappa]*Power2(1/(-psiexp(1,2) + psiexp(2,1)))*psiexp(j - kappa,2)*psiexp(j - k + kappa - nu - 2*Max(0,-j + kappa),1)*(Max(0,-j + kappa)*(-psiexp(1,2) + psiexp(2,1))*psiexp(2*(-j + k),1)*psiexp(1 + Max(0,-j + kappa),2) + (1 - n + q)*psiexp(2 - 2*j + 2*k,1)*psiexp(1 + Max(0,-j + kappa),2) + (-1 - j + k - n + q)*psiexp(1 - j + k,2)*psiexp(2*(1 + Max(0,-j + kappa)),1) + (n - q)*psiexp(2*(-j + k),1)*psiexp(2 + Max(0,-j + kappa),2) + (j - k + n - q)*psiexp(-j + k,2)*psiexp(4 + 2*Max(0,-j + kappa),1));
					}
					sumInuP1 = sumIkappa;
				}
				sumInu += sumInuP1*b[nu];
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += MMjkP1/(n - q);
	}
}

for (int j = 1; j <= -2 + q; j++)
{
	for (int k = Max(3,2 + j); k <= q; k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = Max(0,2 + j - k); nu <= -1 + j; nu++)
			{
				scalar sumInuP1;
				{
					scalar sumIkappa = 0;
					for (int kappa = 0; kappa <= -2 - j + k + nu; kappa++)
					{
						sumIkappa += -(vEG2G2G2G2c*(j - k + 2*n - 2*q - Max(0,kappa - nu))*(1 + j - k + Max(0,kappa - nu))*b[kappa])/2. - (vEG2G2G2G2cAWADAH*(-1 - j + k - Max(0,kappa - nu))*(-j + k - 2*n + 2*q + Max(0,kappa - nu))*b[kappa]*psiexp(k - kappa,1))/2. - (vEG2G2G2G2cAD*(j - k + 2*n - 2*q - Max(0,kappa - nu))*(1 + j - k + Max(0,kappa - nu))*b[kappa]*psiexp(j - nu,1))/2. - (vEG2G2G2G2cAWAH*(j - k + 2*n - 2*q - Max(0,kappa - nu))*(1 + j - k + Max(0,kappa - nu))*b[kappa]*psiexp(-j + k - kappa + nu,1))/2. - (vEG2G2G2G2cAWAHBD*(j - k + 2*n - 2*q - Max(0,kappa - nu))*(1 + j - k + Max(0,kappa - nu))*b[kappa]*psiexp(j - nu,2)*psiexp(-j + k - kappa + nu,1))/2. - vEG2G2G2G2cAH*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-Max(0,kappa - nu),1)*((n - q)*psiexp(-j + k,1) + Max(0,kappa - nu)*(-1 + psiexp(1,1))*psiexp(-j + k,1) + (1 - n + q)*psiexp(1 - j + k,1) + (-1 - j + k - n + q)*psiexp(1 + Max(0,kappa - nu),1) + (j - k + n - q)*psiexp(2 + Max(0,kappa - nu),1)) + vEG2G2G2G2cAWAD*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(j - kappa,1)*((-j + k - n + q)*psiexp(-j + k,1) + (1 + j - k + n - q)*psiexp(1 - j + k,1) + (-1 + n - q)*psiexp(1 + Max(0,kappa - nu),1) + Max(0,kappa - nu)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,kappa - nu),1) + (-n + q)*psiexp(2 + Max(0,kappa - nu),1)) + vEG2G2G2G2cAW*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-kappa + nu,1)*((-j + k - n + q)*psiexp(-j + k,1) + (1 + j - k + n - q)*psiexp(1 - j + k,1) + (-1 + n - q)*psiexp(1 + Max(0,kappa - nu),1) + Max(0,kappa - nu)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,kappa - nu),1) + (-n + q)*psiexp(2 + Max(0,kappa - nu),1)) + vEG2G2G2G2cAWBD*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(j - nu,2)*psiexp(-kappa + nu,1)*((-j + k - n + q)*psiexp(-j + k,1) + (1 + j - k + n - q)*psiexp(1 - j + k,1) + (-1 + n - q)*psiexp(1 + Max(0,kappa - nu),1) + Max(0,kappa - nu)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,kappa - nu),1) + (-n + q)*psiexp(2 + Max(0,kappa - nu),1)) + vEG2G2G2G2cADAH*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(j - nu - Max(0,kappa - nu),1)*((-n + q)*psiexp(-j + k,1) - Max(0,kappa - nu)*(-1 + psiexp(1,1))*psiexp(-j + k,1) + (-1 + n - q)*psiexp(1 - j + k,1) + (1 + j - k + n - q)*psiexp(1 + Max(0,kappa - nu),1) + (-j + k - n + q)*psiexp(2 + Max(0,kappa - nu),1)) + vEG2G2G2G2cADBH*b[kappa]*Power2(1/(-1 + psiexp(1,2)))*psiexp(j - nu,1)*psiexp(-Max(0,kappa - nu),2)*((-n + q)*psiexp(-j + k,2) - Max(0,kappa - nu)*(-1 + psiexp(1,2))*psiexp(-j + k,2) + (-1 + n - q)*psiexp(1 - j + k,2) + (1 + j - k + n - q)*psiexp(1 + Max(0,kappa - nu),2) + (-j + k - n + q)*psiexp(2 + Max(0,kappa - nu),2)) - vEG2G2G2G2cAWADBH*b[kappa]*Power2(1/(psiexp(1,1) - psiexp(1,2)))*psiexp(j - kappa,1)*psiexp(-Max(0,kappa - nu),2)*(-(Max(0,kappa - nu)*(psiexp(1,1) - psiexp(1,2))*psiexp(-j + k,2)*psiexp(1 + Max(0,kappa - nu),1)) + (1 - n + q)*psiexp(1 - j + k,2)*psiexp(1 + Max(0,kappa - nu),1) + (-1 - j + k - n + q)*psiexp(1 - j + k,1)*psiexp(1 + Max(0,kappa - nu),2) + (n - q)*psiexp(-j + k,2)*psiexp(2 + Max(0,kappa - nu),1) + (j - k + n - q)*psiexp(-j + k,1)*psiexp(2 + Max(0,kappa - nu),2)) + vEG2G2G2G2cAWBDBH*b[kappa]*Power2(1/(psiexp(1,1) - psiexp(1,2)))*psiexp(-kappa + nu,1)*psiexp(j - nu - Max(0,kappa - nu),2)*(Max(0,kappa - nu)*(psiexp(1,1) - psiexp(1,2))*psiexp(-j + k,2)*psiexp(1 + Max(0,kappa - nu),1) + (-1 + n - q)*psiexp(1 - j + k,2)*psiexp(1 + Max(0,kappa - nu),1) + (1 + j - k + n - q)*psiexp(1 - j + k,1)*psiexp(1 + Max(0,kappa - nu),2) + (-n + q)*psiexp(-j + k,2)*psiexp(2 + Max(0,kappa - nu),1) + (-j + k - n + q)*psiexp(-j + k,1)*psiexp(2 + Max(0,kappa - nu),2)) - vEG2G2G2G2cAWBDCH*b[kappa]*Power2(1/(psiexp(1,1) - psiexp(1,3)))*psiexp(j - nu,2)*psiexp(-kappa + nu,1)*psiexp(-Max(0,kappa - nu),3)*(-(Max(0,kappa - nu)*(psiexp(1,1) - psiexp(1,3))*psiexp(-j + k,3)*psiexp(1 + Max(0,kappa - nu),1)) + (1 - n + q)*psiexp(1 - j + k,3)*psiexp(1 + Max(0,kappa - nu),1) + (-1 - j + k - n + q)*psiexp(1 - j + k,1)*psiexp(1 + Max(0,kappa - nu),3) + (n - q)*psiexp(-j + k,3)*psiexp(2 + Max(0,kappa - nu),1) + (j - k + n - q)*psiexp(-j + k,1)*psiexp(2 + Max(0,kappa - nu),3));
					}
					sumInuP1 = sumIkappa;
				}
				sumInu += sumInuP1*b[nu];
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += MMjkP1/(n - q);
	}
}

for (int j = 1; j <= -2 + q; j++)
{
	for (int k = Max(1,2 + j); k <= q; k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = 0; nu <= -1 + j; nu++)
			{
				sumInu += -((1 + j - k)*(j - k + 2*n - 2*q)*vEG2G2G4c*b[nu])/2. + vEG2G2G4cAD*b[nu]*Power2(1/(-1 + psiexp(1,1)))*(-n + q + (-1 + n - q)*psiexp(1,1) + (1 + j - k + n - q)*psiexp(1 + j - k,1) + (-j + k - n + q)*psiexp(2 + j - k,1))*psiexp(-j + k,1) - ((1 + j - k)*(j - k + 2*n - 2*q)*vEG2G2G4cAW*b[nu]*psiexp(j - nu,1))/2. + vEG2G2G4cAWBD*b[nu]*Power2(1/(-1 + psiexp(1,2)))*(-n + q + (-1 + n - q)*psiexp(1,2) + (1 + j - k + n - q)*psiexp(1 + j - k,2) + (-j + k - n + q)*psiexp(2 + j - k,2))*psiexp(-j + k,2)*psiexp(j - nu,1) + vEG2G2G4cAWAD*b[nu]*Power2(1/(-1 + psiexp(1,1)))*(-n + q + (-1 + n - q)*psiexp(1,1) + (1 + j - k + n - q)*psiexp(1 + j - k,1) + (-j + k - n + q)*psiexp(2 + j - k,1))*psiexp(k - nu,1);
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += (MMjkP1*b[k])/(n - q);
	}
}

for (int j = 1; j <= -2 + q; j++)
{
	for (int k = Max(1,2 + j); k <= q; k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = 0; nu <= -1 + j; nu++)
			{
				scalar sumInuP1;
				{
					scalar sumIkappa = 0;
					for (int kappa = Max(0,1 + j); kappa <= -1 + k; kappa++)
					{
						sumIkappa += (j - kappa + n - q)*vEG2G4G2c*b[kappa] + (j - kappa + n - q)*vEG2G4G2cAH*b[kappa]*psiexp(k - kappa,1) + (j - kappa + n - q)*vEG2G4G2cBH*b[kappa]*psiexp(k - kappa,2) + (j - kappa + n - q)*vEG2G4G2cAW*b[kappa]*psiexp(j - nu,1) + (j - kappa + n - q)*vEG2G4G2cAWBH*b[kappa]*psiexp(k - kappa,2)*psiexp(j - nu,1) + (j - kappa + n - q)*vEG2G4G2cAWCH*b[kappa]*psiexp(k - kappa,3)*psiexp(j - nu,1) + (j - kappa + n - q)*vEG2G4G2cAWAH*b[kappa]*psiexp(j + k - kappa - nu,1);
					}
					sumInuP1 = sumIkappa;
				}
				sumInu += sumInuP1*b[nu];
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += MMjkP1/(n - q);
	}
}

for (int j = 1; j <= -2 + q; j++)
{
	for (int k = Max(2,2 + j); k <= q; k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = Max(0,1 + j - k); nu <= -1 + j; nu++)
			{
				scalar sumInuP1;
				{
					scalar sumIkappa = 0;
					for (int kappa = Max(0,1 + nu); kappa <= Min(-1 - j + k + nu,-1 + n + nu - q); kappa++)
					{
						sumIkappa += (-kappa + n + nu - q)*vEG4G2G2c*b[kappa] - (kappa - n - nu + q)*vEG4G2G2cADAH*b[kappa]*psiexp(k - kappa,1) - (kappa - n - nu + q)*vEG4G2G2cBDBH*b[kappa]*psiexp(k - kappa,2) + (-kappa + n + nu - q)*vEG4G2G2cAD*b[kappa]*psiexp(j - nu,1) + (-kappa + n + nu - q)*vEG4G2G2cBD*b[kappa]*psiexp(j - nu,2) - (kappa - n - nu + q)*vEG4G2G2cAH*b[kappa]*psiexp(-j + k - kappa + nu,1) + (-kappa + n + nu - q)*vEG4G2G2cAHBD*b[kappa]*psiexp(j - nu,2)*psiexp(-j + k - kappa + nu,1) + (-kappa + n + nu - q)*vEG4G2G2cADBH*b[kappa]*psiexp(j - nu,1)*psiexp(-j + k - kappa + nu,2) + (-kappa + n + nu - q)*vEG4G2G2cBDCH*b[kappa]*psiexp(j - nu,2)*psiexp(-j + k - kappa + nu,3);
					}
					sumInuP1 = sumIkappa;
				}
				sumInu += sumInuP1*b[nu];
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += MMjkP1/(n - q);
	}
}

for (int j = 2; j <= -3 + q; j++)
{
	for (int k = 2; k <= -3 + q; k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = Max(0,1 + j + k - q); nu <= Min(-2 + j,-2 + k); nu++)
			{
				scalar sumInuP1;
				{
					scalar sumIkappa = 0;
					for (int kappa = 1 + j + k - nu; kappa <= q; kappa++)
					{
						sumIkappa += -(vEG2G2G2G2c*(-1 + k - nu - Max(0,-j + k))*(k - 2*n - nu + 2*q + Max(0,-j + k))*b[kappa])/2. - (vEG2G2G2G2cADAH*(-1 + k - nu - Max(0,-j + k))*(k - 2*n - nu + 2*q + Max(0,-j + k))*b[kappa]*psiexp(-k + kappa,1))/2. - (vEG2G2G2G2cAWAD*(-k + 2*n + nu - 2*q - Max(0,-j + k))*(1 - k + nu + Max(0,-j + k))*b[kappa]*psiexp(j - nu,1))/2. + vEG2G2G2G2cAD*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(j - k,1)*((k - n - nu + q)*psiexp(k - nu,1) + (1 - k + n + nu - q)*psiexp(1 + k - nu,1) + (-1 + n - q)*psiexp(1 + Max(0,-j + k),1) + Max(0,-j + k)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,-j + k),1) + (-n + q)*psiexp(2 + Max(0,-j + k),1)) + vEG2G2G2G2cAWADAH*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-k + kappa - Max(0,-j + k),1)*((-n + q)*psiexp(k - nu,1) - Max(0,-j + k)*(-1 + psiexp(1,1))*psiexp(k - nu,1) + (-1 + n - q)*psiexp(1 + k - nu,1) + (1 - k + n + nu - q)*psiexp(1 + Max(0,-j + k),1) + (k - n - nu + q)*psiexp(2 + Max(0,-j + k),1)) + vEG2G2G2G2cAH*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-j - k + kappa + nu - Max(0,-j + k),1)*((-n + q)*psiexp(k - nu,1) - Max(0,-j + k)*(-1 + psiexp(1,1))*psiexp(k - nu,1) + (-1 + n - q)*psiexp(1 + k - nu,1) + (1 - k + n + nu - q)*psiexp(1 + Max(0,-j + k),1) + (k - n - nu + q)*psiexp(2 + Max(0,-j + k),1)) + vEG2G2G2G2cAW*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-Max(0,-j + k),1)*((-n + q)*psiexp(k - nu,1) - Max(0,-j + k)*(-1 + psiexp(1,1))*psiexp(k - nu,1) + (-1 + n - q)*psiexp(1 + k - nu,1) + (1 - k + n + nu - q)*psiexp(1 + Max(0,-j + k),1) + (k - n - nu + q)*psiexp(2 + Max(0,-j + k),1)) + vEG2G2G2G2cAWBDBH*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-k + kappa,2)*psiexp(-Max(0,-j + k),1)*((-n + q)*psiexp(k - nu,1) - Max(0,-j + k)*(-1 + psiexp(1,1))*psiexp(k - nu,1) + (-1 + n - q)*psiexp(1 + k - nu,1) + (1 - k + n + nu - q)*psiexp(1 + Max(0,-j + k),1) + (k - n - nu + q)*psiexp(2 + Max(0,-j + k),1)) - vEG2G2G2G2cAWADBH*b[kappa]*Power2(1/(-1 + psiexp(1,2)))*psiexp(j - nu,1)*psiexp(-j - k + kappa + nu - Max(0,-j + k),2)*((n - q)*psiexp(k - nu,2) + Max(0,-j + k)*(-1 + psiexp(1,2))*psiexp(k - nu,2) + (1 - n + q)*psiexp(1 + k - nu,2) + (-1 + k - n - nu + q)*psiexp(1 + Max(0,-j + k),2) + (-k + n + nu - q)*psiexp(2 + Max(0,-j + k),2)) + vEG2G2G2G2cAWBD*b[kappa]*Power2(1/(psiexp(1,1) - psiexp(1,2)))*psiexp(j - k,2)*psiexp(-Max(0,-j + k),1)*((1 - k + n + nu - q)*psiexp(1 + k - nu,2)*psiexp(1 + Max(0,-j + k),1) + Max(0,-j + k)*(-psiexp(1,1) + psiexp(1,2))*psiexp(k - nu,1)*psiexp(1 + Max(0,-j + k),2) + (-1 + n - q)*psiexp(1 + k - nu,1)*psiexp(1 + Max(0,-j + k),2) + (k - n - nu + q)*psiexp(k - nu,2)*psiexp(2 + Max(0,-j + k),1) + (-n + q)*psiexp(k - nu,1)*psiexp(2 + Max(0,-j + k),2)) + vEG2G2G2G2cADBH*b[kappa]*Power2(1/(psiexp(1,1) - psiexp(1,2)))*psiexp(j - k,1)*psiexp(-j - k + kappa + nu - Max(0,-j + k),2)*(Max(0,-j + k)*(psiexp(1,1) - psiexp(1,2))*psiexp(k - nu,2)*psiexp(1 + Max(0,-j + k),1) + (-1 + n - q)*psiexp(1 + k - nu,2)*psiexp(1 + Max(0,-j + k),1) + (1 - k + n + nu - q)*psiexp(1 + k - nu,1)*psiexp(1 + Max(0,-j + k),2) + (-n + q)*psiexp(k - nu,2)*psiexp(2 + Max(0,-j + k),1) + (k - n - nu + q)*psiexp(k - nu,1)*psiexp(2 + Max(0,-j + k),2)) + vEG2G2G2G2cAWAH*b[kappa]*Power2(1/(-1 + psiexp(2,1)))*psiexp(-j - k + kappa + nu - 2*Max(0,-j + k),1)*((-n + q)*psiexp(2*(k - nu),1) - Max(0,-j + k)*(-1 + psiexp(2,1))*psiexp(2*(k - nu),1) + (-1 + n - q)*psiexp(2*(1 + k - nu),1) + (1 - k + n + nu - q)*psiexp(2*(1 + Max(0,-j + k)),1) + (k - n - nu + q)*psiexp(2*(2 + Max(0,-j + k)),1)) + vEG2G2G2G2cAWBDCH*b[kappa]*Power2(1/(psiexp(1,2) - psiexp(1,1)*psiexp(1,3)))*psiexp(j - k - nu,2)*psiexp(-j - k + kappa - Max(0,-j + k),3)*psiexp(-nu - Max(0,-j + k),1)*(Max(0,-j + k)*(psiexp(1,2) - psiexp(1,1)*psiexp(1,3))*psiexp(k,1)*psiexp(k,3)*psiexp(1 + nu + Max(0,-j + k),2) + (-1 + n - q)*psiexp(1 + k,1)*psiexp(1 + k,3)*psiexp(1 + nu + Max(0,-j + k),2) + (1 - k + n + nu - q)*psiexp(1 + k,2)*psiexp(1 + nu + Max(0,-j + k),1)*psiexp(1 + nu + Max(0,-j + k),3) + (-n + q)*psiexp(k,1)*psiexp(k,3)*psiexp(2 + nu + Max(0,-j + k),2) + (k - n - nu + q)*psiexp(k,2)*psiexp(2 + nu + Max(0,-j + k),1)*psiexp(2 + nu + Max(0,-j + k),3)) + vEG2G2G2G2cAWAHBD*b[kappa]*Power2(1/(-psiexp(1,2) + psiexp(2,1)))*psiexp(j - k,2)*psiexp(-j - k + kappa + nu - 2*Max(0,-j + k),1)*(Max(0,-j + k)*(psiexp(1,2) - psiexp(2,1))*psiexp(2*k - 2*nu,1)*psiexp(1 + Max(0,-j + k),2) + (-1 + n - q)*psiexp(2*(1 + k - nu),1)*psiexp(1 + Max(0,-j + k),2) + (1 - k + n + nu - q)*psiexp(1 + k - nu,2)*psiexp(2*(1 + Max(0,-j + k)),1) + (-n + q)*psiexp(2*(k - nu),1)*psiexp(2 + Max(0,-j + k),2) + (k - n - nu + q)*psiexp(k - nu,2)*psiexp(4 + 2*Max(0,-j + k),1));
					}
					sumInuP1 = sumIkappa;
				}
				sumInu += sumInuP1*b[nu];
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += MMjkP1/(n - q);
	}
}

for (int j = 2; j <= -2 + q; j++)
{
	for (int k = 2; k <= -2 + q; k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = 0; nu <= Min(-2 + j,-2 + k); nu++)
			{
				scalar sumInuP1;
				{
					scalar sumIkappa = 0;
					for (int kappa = Max(2 + j,2 + k); kappa <= Min(j + k - nu,q); kappa++)
					{
						sumIkappa += -(vEG2G2G2G2c*(j - kappa + 2*n - 2*q - Max(0,-j + k))*(1 + j - kappa + Max(0,-j + k))*b[kappa])/2. - (vEG2G2G2G2cADAH*(j - kappa + 2*n - 2*q - Max(0,-j + k))*(1 + j - kappa + Max(0,-j + k))*b[kappa]*psiexp(-k + kappa,1))/2. - (vEG2G2G2G2cAWAD*(j - kappa + 2*n - 2*q - Max(0,-j + k))*(1 + j - kappa + Max(0,-j + k))*b[kappa]*psiexp(j - nu,1))/2. - vEG2G2G2G2cAW*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(j + k - kappa - nu - Max(0,-j + k),1)*((n - q)*psiexp(-j + kappa,1) + Max(0,-j + k)*(-1 + psiexp(1,1))*psiexp(-j + kappa,1) + (1 - n + q)*psiexp(1 - j + kappa,1) + (-1 - j + kappa - n + q)*psiexp(1 + Max(0,-j + k),1) + (j - kappa + n - q)*psiexp(2 + Max(0,-j + k),1)) - vEG2G2G2G2cAH*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-Max(0,-j + k),1)*((n - q)*psiexp(-j + kappa,1) + Max(0,-j + k)*(-1 + psiexp(1,1))*psiexp(-j + kappa,1) + (1 - n + q)*psiexp(1 - j + kappa,1) + (-1 - j + kappa - n + q)*psiexp(1 + Max(0,-j + k),1) + (j - kappa + n - q)*psiexp(2 + Max(0,-j + k),1)) + vEG2G2G2G2cAD*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(j - k,1)*((-j + kappa - n + q)*psiexp(-j + kappa,1) + (1 + j - kappa + n - q)*psiexp(1 - j + kappa,1) + (-1 + n - q)*psiexp(1 + Max(0,-j + k),1) + Max(0,-j + k)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,-j + k),1) + (-n + q)*psiexp(2 + Max(0,-j + k),1)) + vEG2G2G2G2cAWADAH*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(j - nu - Max(0,-j + k),1)*((-n + q)*psiexp(-j + kappa,1) - Max(0,-j + k)*(-1 + psiexp(1,1))*psiexp(-j + kappa,1) + (-1 + n - q)*psiexp(1 - j + kappa,1) + (1 + j - kappa + n - q)*psiexp(1 + Max(0,-j + k),1) + (-j + kappa - n + q)*psiexp(2 + Max(0,-j + k),1)) + vEG2G2G2G2cAWBDBH*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-k + kappa,2)*psiexp(j + k - kappa - nu - Max(0,-j + k),1)*((-n + q)*psiexp(-j + kappa,1) - Max(0,-j + k)*(-1 + psiexp(1,1))*psiexp(-j + kappa,1) + (-1 + n - q)*psiexp(1 - j + kappa,1) + (1 + j - kappa + n - q)*psiexp(1 + Max(0,-j + k),1) + (-j + kappa - n + q)*psiexp(2 + Max(0,-j + k),1)) + vEG2G2G2G2cAWADBH*b[kappa]*Power2(1/(-1 + psiexp(1,2)))*psiexp(j - nu,1)*psiexp(-Max(0,-j + k),2)*((-n + q)*psiexp(-j + kappa,2) - Max(0,-j + k)*(-1 + psiexp(1,2))*psiexp(-j + kappa,2) + (-1 + n - q)*psiexp(1 - j + kappa,2) + (1 + j - kappa + n - q)*psiexp(1 + Max(0,-j + k),2) + (-j + kappa - n + q)*psiexp(2 + Max(0,-j + k),2)) - vEG2G2G2G2cADBH*b[kappa]*Power2(1/(psiexp(1,1) - psiexp(1,2)))*psiexp(j - k,1)*psiexp(-Max(0,-j + k),2)*(-(Max(0,-j + k)*(psiexp(1,1) - psiexp(1,2))*psiexp(-j + kappa,2)*psiexp(1 + Max(0,-j + k),1)) + (1 - n + q)*psiexp(1 - j + kappa,2)*psiexp(1 + Max(0,-j + k),1) + (-1 - j + kappa - n + q)*psiexp(1 - j + kappa,1)*psiexp(1 + Max(0,-j + k),2) + (n - q)*psiexp(-j + kappa,2)*psiexp(2 + Max(0,-j + k),1) + (j - kappa + n - q)*psiexp(-j + kappa,1)*psiexp(2 + Max(0,-j + k),2)) + vEG2G2G2G2cAWBD*b[kappa]*Power2(1/(psiexp(1,1) - psiexp(1,2)))*psiexp(j - k,2)*psiexp(j + k - kappa - nu - Max(0,-j + k),1)*((1 + j - kappa + n - q)*psiexp(1 - j + kappa,2)*psiexp(1 + Max(0,-j + k),1) + Max(0,-j + k)*(-psiexp(1,1) + psiexp(1,2))*psiexp(-j + kappa,1)*psiexp(1 + Max(0,-j + k),2) + (-1 + n - q)*psiexp(1 - j + kappa,1)*psiexp(1 + Max(0,-j + k),2) + (-j + kappa - n + q)*psiexp(-j + kappa,2)*psiexp(2 + Max(0,-j + k),1) + (-n + q)*psiexp(-j + kappa,1)*psiexp(2 + Max(0,-j + k),2)) + vEG2G2G2G2cAWAH*b[kappa]*Power2(1/(-1 + psiexp(2,1)))*psiexp(j + k - kappa - nu - 2*Max(0,-j + k),1)*((-n + q)*psiexp(2*(-j + kappa),1) - Max(0,-j + k)*(-1 + psiexp(2,1))*psiexp(2*(-j + kappa),1) + (-1 + n - q)*psiexp(2*(1 - j + kappa),1) + (1 + j - kappa + n - q)*psiexp(2*(1 + Max(0,-j + k)),1) + (-j + kappa - n + q)*psiexp(2*(2 + Max(0,-j + k)),1)) - vEG2G2G2G2cAWBDCH*b[kappa]*Power2(1/(psiexp(1,2) - psiexp(1,1)*psiexp(1,3)))*psiexp(-k,2)*psiexp(-j - Max(0,-j + k),3)*psiexp(k - kappa - nu - Max(0,-j + k),1)*(Max(0,-j + k)*(-psiexp(1,2) + psiexp(1,1)*psiexp(1,3))*psiexp(kappa,1)*psiexp(kappa,3)*psiexp(1 + j + Max(0,-j + k),2) + (1 - n + q)*psiexp(1 + kappa,1)*psiexp(1 + kappa,3)*psiexp(1 + j + Max(0,-j + k),2) + (-1 - j + kappa - n + q)*psiexp(1 + kappa,2)*psiexp(1 + j + Max(0,-j + k),1)*psiexp(1 + j + Max(0,-j + k),3) + (n - q)*psiexp(kappa,1)*psiexp(kappa,3)*psiexp(2 + j + Max(0,-j + k),2) + (j - kappa + n - q)*psiexp(kappa,2)*psiexp(2 + j + Max(0,-j + k),1)*psiexp(2 + j + Max(0,-j + k),3)) - vEG2G2G2G2cAWAHBD*b[kappa]*Power2(1/(-psiexp(1,2) + psiexp(2,1)))*psiexp(j - k,2)*psiexp(j + k - kappa - nu - 2*Max(0,-j + k),1)*(Max(0,-j + k)*(-psiexp(1,2) + psiexp(2,1))*psiexp(2*(-j + kappa),1)*psiexp(1 + Max(0,-j + k),2) + (1 - n + q)*psiexp(2 - 2*j + 2*kappa,1)*psiexp(1 + Max(0,-j + k),2) + (-1 - j + kappa - n + q)*psiexp(1 - j + kappa,2)*psiexp(2*(1 + Max(0,-j + k)),1) + (n - q)*psiexp(2*(-j + kappa),1)*psiexp(2 + Max(0,-j + k),2) + (j - kappa + n - q)*psiexp(-j + kappa,2)*psiexp(4 + 2*Max(0,-j + k),1));
					}
					sumInuP1 = sumIkappa;
				}
				sumInu += sumInuP1*b[nu];
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += MMjkP1/(n - q);
	}
}

for (int j = 3; j <= q; j++)
{
	for (int k = 2; k <= -1 + q; k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = 0; nu <= Min(-3 + j,-2 + k); nu++)
			{
				scalar sumInuP1;
				{
					scalar sumIkappa = 0;
					for (int kappa = 1 + k; kappa <= Min(-2 + j + k - nu,q); kappa++)
					{
						sumIkappa += -(vEG2G2G2G2c*(-1 + k - nu - Max(0,-j + kappa))*(k - 2*n - nu + 2*q + Max(0,-j + kappa))*b[kappa])/2. - (vEG2G2G2G2cAD*(-1 + k - nu - Max(0,-j + kappa))*(k - 2*n - nu + 2*q + Max(0,-j + kappa))*b[kappa]*psiexp(-k + kappa,1))/2. - (vEG2G2G2G2cAWADAH*(-k + 2*n + nu - 2*q - Max(0,-j + kappa))*(1 - k + nu + Max(0,-j + kappa))*b[kappa]*psiexp(j - nu,1))/2. - (vEG2G2G2G2cAWAH*(-1 + k - nu - Max(0,-j + kappa))*(k - 2*n - nu + 2*q + Max(0,-j + kappa))*b[kappa]*psiexp(j + k - kappa - nu,1))/2. - (vEG2G2G2G2cAWAHBD*(-1 + k - nu - Max(0,-j + kappa))*(k - 2*n - nu + 2*q + Max(0,-j + kappa))*b[kappa]*psiexp(-k + kappa,2)*psiexp(j + k - kappa - nu,1))/2. + vEG2G2G2G2cADAH*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(j - k,1)*((k - n - nu + q)*psiexp(k - nu,1) + (1 - k + n + nu - q)*psiexp(1 + k - nu,1) + (-1 + n - q)*psiexp(1 + Max(0,-j + kappa),1) + Max(0,-j + kappa)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,-j + kappa),1) + (-n + q)*psiexp(2 + Max(0,-j + kappa),1)) + vEG2G2G2G2cAH*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(j - kappa,1)*((k - n - nu + q)*psiexp(k - nu,1) + (1 - k + n + nu - q)*psiexp(1 + k - nu,1) + (-1 + n - q)*psiexp(1 + Max(0,-j + kappa),1) + Max(0,-j + kappa)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,-j + kappa),1) + (-n + q)*psiexp(2 + Max(0,-j + kappa),1)) + vEG2G2G2G2cAWAD*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-k + kappa - Max(0,-j + kappa),1)*((-n + q)*psiexp(k - nu,1) - Max(0,-j + kappa)*(-1 + psiexp(1,1))*psiexp(k - nu,1) + (-1 + n - q)*psiexp(1 + k - nu,1) + (1 - k + n + nu - q)*psiexp(1 + Max(0,-j + kappa),1) + (k - n - nu + q)*psiexp(2 + Max(0,-j + kappa),1)) + vEG2G2G2G2cAW*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-Max(0,-j + kappa),1)*((-n + q)*psiexp(k - nu,1) - Max(0,-j + kappa)*(-1 + psiexp(1,1))*psiexp(k - nu,1) + (-1 + n - q)*psiexp(1 + k - nu,1) + (1 - k + n + nu - q)*psiexp(1 + Max(0,-j + kappa),1) + (k - n - nu + q)*psiexp(2 + Max(0,-j + kappa),1)) + vEG2G2G2G2cAWBD*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-k + kappa,2)*psiexp(-Max(0,-j + kappa),1)*((-n + q)*psiexp(k - nu,1) - Max(0,-j + kappa)*(-1 + psiexp(1,1))*psiexp(k - nu,1) + (-1 + n - q)*psiexp(1 + k - nu,1) + (1 - k + n + nu - q)*psiexp(1 + Max(0,-j + kappa),1) + (k - n - nu + q)*psiexp(2 + Max(0,-j + kappa),1)) + vEG2G2G2G2cADBH*b[kappa]*Power2(1/(-1 + psiexp(1,2)))*psiexp(j - kappa,2)*psiexp(-k + kappa,1)*((k - n - nu + q)*psiexp(k - nu,2) + (1 - k + n + nu - q)*psiexp(1 + k - nu,2) + (-1 + n - q)*psiexp(1 + Max(0,-j + kappa),2) + Max(0,-j + kappa)*(-1 + psiexp(1,2))*psiexp(1 + Max(0,-j + kappa),2) + (-n + q)*psiexp(2 + Max(0,-j + kappa),2)) + vEG2G2G2G2cAWADBH*b[kappa]*Power2(1/(psiexp(1,1) - psiexp(1,2)))*psiexp(j - kappa,2)*psiexp(-k + kappa - Max(0,-j + kappa),1)*((1 - k + n + nu - q)*psiexp(1 + k - nu,2)*psiexp(1 + Max(0,-j + kappa),1) + Max(0,-j + kappa)*(-psiexp(1,1) + psiexp(1,2))*psiexp(k - nu,1)*psiexp(1 + Max(0,-j + kappa),2) + (-1 + n - q)*psiexp(1 + k - nu,1)*psiexp(1 + Max(0,-j + kappa),2) + (k - n - nu + q)*psiexp(k - nu,2)*psiexp(2 + Max(0,-j + kappa),1) + (-n + q)*psiexp(k - nu,1)*psiexp(2 + Max(0,-j + kappa),2)) + vEG2G2G2G2cAWBDBH*b[kappa]*Power2(1/(psiexp(1,1) - psiexp(1,2)))*psiexp(j - k,2)*psiexp(-Max(0,-j + kappa),1)*((1 - k + n + nu - q)*psiexp(1 + k - nu,2)*psiexp(1 + Max(0,-j + kappa),1) + Max(0,-j + kappa)*(-psiexp(1,1) + psiexp(1,2))*psiexp(k - nu,1)*psiexp(1 + Max(0,-j + kappa),2) + (-1 + n - q)*psiexp(1 + k - nu,1)*psiexp(1 + Max(0,-j + kappa),2) + (k - n - nu + q)*psiexp(k - nu,2)*psiexp(2 + Max(0,-j + kappa),1) + (-n + q)*psiexp(k - nu,1)*psiexp(2 + Max(0,-j + kappa),2)) + vEG2G2G2G2cAWBDCH*b[kappa]*Power2(1/(psiexp(1,1) - psiexp(1,3)))*psiexp(j - kappa,3)*psiexp(-k + kappa,2)*psiexp(-Max(0,-j + kappa),1)*((1 - k + n + nu - q)*psiexp(1 + k - nu,3)*psiexp(1 + Max(0,-j + kappa),1) + Max(0,-j + kappa)*(-psiexp(1,1) + psiexp(1,3))*psiexp(k - nu,1)*psiexp(1 + Max(0,-j + kappa),3) + (-1 + n - q)*psiexp(1 + k - nu,1)*psiexp(1 + Max(0,-j + kappa),3) + (k - n - nu + q)*psiexp(k - nu,3)*psiexp(2 + Max(0,-j + kappa),1) + (-n + q)*psiexp(k - nu,1)*psiexp(2 + Max(0,-j + kappa),3));
					}
					sumInuP1 = sumIkappa;
				}
				sumInu += sumInuP1*b[nu];
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += MMjkP1/(n - q);
	}
}

for (int j = 3; j <= q; j++)
{
	for (int k = 3; k <= q; k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = 0; nu <= Min(-3 + j,-3 + k); nu++)
			{
				scalar sumInuP1;
				{
					scalar sumIkappa = 0;
					for (int kappa = Max3(0,2 + nu,2 - j + k + nu); kappa <= -1 + k; kappa++)
					{
						sumIkappa += -(vEG2G2G2G2c*(-1 + kappa - nu - Max(0,-j + k))*(kappa - 2*n - nu + 2*q + Max(0,-j + k))*b[kappa])/2. - (vEG2G2G2G2cAD*(-1 + kappa - nu - Max(0,-j + k))*(kappa - 2*n - nu + 2*q + Max(0,-j + k))*b[kappa]*psiexp(k - kappa,1))/2. - (vEG2G2G2G2cAWADAH*(-kappa + 2*n + nu - 2*q - Max(0,-j + k))*(1 - kappa + nu + Max(0,-j + k))*b[kappa]*psiexp(j - nu,1))/2. - (vEG2G2G2G2cAWAH*(-1 + kappa - nu - Max(0,-j + k))*(kappa - 2*n - nu + 2*q + Max(0,-j + k))*b[kappa]*psiexp(j - k + kappa - nu,1))/2. - (vEG2G2G2G2cAWAHBD*(-1 + kappa - nu - Max(0,-j + k))*(kappa - 2*n - nu + 2*q + Max(0,-j + k))*b[kappa]*psiexp(k - kappa,2)*psiexp(j - k + kappa - nu,1))/2. + vEG2G2G2G2cAH*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(j - k,1)*((kappa - n - nu + q)*psiexp(kappa - nu,1) + (1 - kappa + n + nu - q)*psiexp(1 + kappa - nu,1) + (-1 + n - q)*psiexp(1 + Max(0,-j + k),1) + Max(0,-j + k)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,-j + k),1) + (-n + q)*psiexp(2 + Max(0,-j + k),1)) + vEG2G2G2G2cADAH*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(j - kappa,1)*((kappa - n - nu + q)*psiexp(kappa - nu,1) + (1 - kappa + n + nu - q)*psiexp(1 + kappa - nu,1) + (-1 + n - q)*psiexp(1 + Max(0,-j + k),1) + Max(0,-j + k)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,-j + k),1) + (-n + q)*psiexp(2 + Max(0,-j + k),1)) + vEG2G2G2G2cAWAD*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(k - kappa - Max(0,-j + k),1)*((-n + q)*psiexp(kappa - nu,1) - Max(0,-j + k)*(-1 + psiexp(1,1))*psiexp(kappa - nu,1) + (-1 + n - q)*psiexp(1 + kappa - nu,1) + (1 - kappa + n + nu - q)*psiexp(1 + Max(0,-j + k),1) + (kappa - n - nu + q)*psiexp(2 + Max(0,-j + k),1)) + vEG2G2G2G2cAW*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-Max(0,-j + k),1)*((-n + q)*psiexp(kappa - nu,1) - Max(0,-j + k)*(-1 + psiexp(1,1))*psiexp(kappa - nu,1) + (-1 + n - q)*psiexp(1 + kappa - nu,1) + (1 - kappa + n + nu - q)*psiexp(1 + Max(0,-j + k),1) + (kappa - n - nu + q)*psiexp(2 + Max(0,-j + k),1)) + vEG2G2G2G2cAWBD*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(k - kappa,2)*psiexp(-Max(0,-j + k),1)*((-n + q)*psiexp(kappa - nu,1) - Max(0,-j + k)*(-1 + psiexp(1,1))*psiexp(kappa - nu,1) + (-1 + n - q)*psiexp(1 + kappa - nu,1) + (1 - kappa + n + nu - q)*psiexp(1 + Max(0,-j + k),1) + (kappa - n - nu + q)*psiexp(2 + Max(0,-j + k),1)) + vEG2G2G2G2cADBH*b[kappa]*Power2(1/(-1 + psiexp(1,2)))*psiexp(j - k,2)*psiexp(k - kappa,1)*((kappa - n - nu + q)*psiexp(kappa - nu,2) + (1 - kappa + n + nu - q)*psiexp(1 + kappa - nu,2) + (-1 + n - q)*psiexp(1 + Max(0,-j + k),2) + Max(0,-j + k)*(-1 + psiexp(1,2))*psiexp(1 + Max(0,-j + k),2) + (-n + q)*psiexp(2 + Max(0,-j + k),2)) + vEG2G2G2G2cAWADBH*b[kappa]*Power2(1/(psiexp(1,1) - psiexp(1,2)))*psiexp(j - k,2)*psiexp(k - kappa - Max(0,-j + k),1)*((1 - kappa + n + nu - q)*psiexp(1 + kappa - nu,2)*psiexp(1 + Max(0,-j + k),1) + Max(0,-j + k)*(-psiexp(1,1) + psiexp(1,2))*psiexp(kappa - nu,1)*psiexp(1 + Max(0,-j + k),2) + (-1 + n - q)*psiexp(1 + kappa - nu,1)*psiexp(1 + Max(0,-j + k),2) + (kappa - n - nu + q)*psiexp(kappa - nu,2)*psiexp(2 + Max(0,-j + k),1) + (-n + q)*psiexp(kappa - nu,1)*psiexp(2 + Max(0,-j + k),2)) + vEG2G2G2G2cAWBDBH*b[kappa]*Power2(1/(psiexp(1,1) - psiexp(1,2)))*psiexp(j - kappa,2)*psiexp(-Max(0,-j + k),1)*((1 - kappa + n + nu - q)*psiexp(1 + kappa - nu,2)*psiexp(1 + Max(0,-j + k),1) + Max(0,-j + k)*(-psiexp(1,1) + psiexp(1,2))*psiexp(kappa - nu,1)*psiexp(1 + Max(0,-j + k),2) + (-1 + n - q)*psiexp(1 + kappa - nu,1)*psiexp(1 + Max(0,-j + k),2) + (kappa - n - nu + q)*psiexp(kappa - nu,2)*psiexp(2 + Max(0,-j + k),1) + (-n + q)*psiexp(kappa - nu,1)*psiexp(2 + Max(0,-j + k),2)) + vEG2G2G2G2cAWBDCH*b[kappa]*Power2(1/(psiexp(1,1) - psiexp(1,3)))*psiexp(j - k,3)*psiexp(k - kappa,2)*psiexp(-Max(0,-j + k),1)*((1 - kappa + n + nu - q)*psiexp(1 + kappa - nu,3)*psiexp(1 + Max(0,-j + k),1) + Max(0,-j + k)*(-psiexp(1,1) + psiexp(1,3))*psiexp(kappa - nu,1)*psiexp(1 + Max(0,-j + k),3) + (-1 + n - q)*psiexp(1 + kappa - nu,1)*psiexp(1 + Max(0,-j + k),3) + (kappa - n - nu + q)*psiexp(kappa - nu,3)*psiexp(2 + Max(0,-j + k),1) + (-n + q)*psiexp(kappa - nu,1)*psiexp(2 + Max(0,-j + k),3));
					}
					sumInuP1 = sumIkappa;
				}
				sumInu += sumInuP1*b[nu];
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += MMjkP1/(n - q);
	}
}

for (int j = 2; j <= q; j++)
{
	for (int k = Ceiling(Max(2,(2 + j)/2.)); k <= q; k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = Max(0,j - k); nu <= Min(-2 + j,-2 + k); nu++)
			{
				scalar sumInuP1;
				{
					scalar sumIkappa = 0;
					for (int kappa = 0; kappa <= -j + k + nu; kappa++)
					{
						sumIkappa += -(vEG2G2G2G2c*(-1 + k - nu - Max(0,-j + k))*(k - 2*n - nu + 2*q + Max(0,-j + k))*b[kappa])/2. - (vEG2G2G2G2cAWAD*(-1 + k - nu - Max(0,-j + k))*(k - 2*n - nu + 2*q + Max(0,-j + k))*b[kappa]*psiexp(k - kappa,1))/2. - (vEG2G2G2G2cADAH*(-k + 2*n + nu - 2*q - Max(0,-j + k))*(1 - k + nu + Max(0,-j + k))*b[kappa]*psiexp(j - nu,1))/2. + vEG2G2G2G2cAH*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(j - k,1)*((k - n - nu + q)*psiexp(k - nu,1) + (1 - k + n + nu - q)*psiexp(1 + k - nu,1) + (-1 + n - q)*psiexp(1 + Max(0,-j + k),1) + Max(0,-j + k)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,-j + k),1) + (-n + q)*psiexp(2 + Max(0,-j + k),1)) + vEG2G2G2G2cAWADAH*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(j - kappa,1)*((k - n - nu + q)*psiexp(k - nu,1) + (1 - k + n + nu - q)*psiexp(1 + k - nu,1) + (-1 + n - q)*psiexp(1 + Max(0,-j + k),1) + Max(0,-j + k)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,-j + k),1) + (-n + q)*psiexp(2 + Max(0,-j + k),1)) + vEG2G2G2G2cAW*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-kappa + nu,1)*((k - n - nu + q)*psiexp(k - nu,1) + (1 - k + n + nu - q)*psiexp(1 + k - nu,1) + (-1 + n - q)*psiexp(1 + Max(0,-j + k),1) + Max(0,-j + k)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,-j + k),1) + (-n + q)*psiexp(2 + Max(0,-j + k),1)) + vEG2G2G2G2cAWBDBH*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(j - nu,2)*psiexp(-kappa + nu,1)*((k - n - nu + q)*psiexp(k - nu,1) + (1 - k + n + nu - q)*psiexp(1 + k - nu,1) + (-1 + n - q)*psiexp(1 + Max(0,-j + k),1) + Max(0,-j + k)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,-j + k),1) + (-n + q)*psiexp(2 + Max(0,-j + k),1)) + vEG2G2G2G2cAD*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-Max(0,-j + k),1)*((-n + q)*psiexp(k - nu,1) - Max(0,-j + k)*(-1 + psiexp(1,1))*psiexp(k - nu,1) + (-1 + n - q)*psiexp(1 + k - nu,1) + (1 - k + n + nu - q)*psiexp(1 + Max(0,-j + k),1) + (k - n - nu + q)*psiexp(2 + Max(0,-j + k),1)) + vEG2G2G2G2cAWADBH*b[kappa]*Power2(1/(-1 + psiexp(1,2)))*psiexp(j - k,2)*psiexp(k - kappa,1)*((k - n - nu + q)*psiexp(k - nu,2) + (1 - k + n + nu - q)*psiexp(1 + k - nu,2) + (-1 + n - q)*psiexp(1 + Max(0,-j + k),2) + Max(0,-j + k)*(-1 + psiexp(1,2))*psiexp(1 + Max(0,-j + k),2) + (-n + q)*psiexp(2 + Max(0,-j + k),2)) - vEG2G2G2G2cAWBD*b[kappa]*Power2(1/(psiexp(1,1) - psiexp(1,2)))*psiexp(-kappa + nu,1)*psiexp(-Max(0,-j + k),2)*(-(Max(0,-j + k)*(psiexp(1,1) - psiexp(1,2))*psiexp(k - nu,2)*psiexp(1 + Max(0,-j + k),1)) + (1 - n + q)*psiexp(1 + k - nu,2)*psiexp(1 + Max(0,-j + k),1) + (-1 + k - n - nu + q)*psiexp(1 + k - nu,1)*psiexp(1 + Max(0,-j + k),2) + (n - q)*psiexp(k - nu,2)*psiexp(2 + Max(0,-j + k),1) + (-k + n + nu - q)*psiexp(k - nu,1)*psiexp(2 + Max(0,-j + k),2)) + vEG2G2G2G2cADBH*b[kappa]*Power2(1/(psiexp(1,1) - psiexp(1,2)))*psiexp(j - k,2)*psiexp(-Max(0,-j + k),1)*((1 - k + n + nu - q)*psiexp(1 + k - nu,2)*psiexp(1 + Max(0,-j + k),1) + Max(0,-j + k)*(-psiexp(1,1) + psiexp(1,2))*psiexp(k - nu,1)*psiexp(1 + Max(0,-j + k),2) + (-1 + n - q)*psiexp(1 + k - nu,1)*psiexp(1 + Max(0,-j + k),2) + (k - n - nu + q)*psiexp(k - nu,2)*psiexp(2 + Max(0,-j + k),1) + (-n + q)*psiexp(k - nu,1)*psiexp(2 + Max(0,-j + k),2)) + vEG2G2G2G2cAWAH*b[kappa]*Power2(1/(-1 + psiexp(2,1)))*psiexp(j - k - kappa + nu,1)*((k - n - nu + q)*psiexp(2*(k - nu),1) + (1 - k + n + nu - q)*psiexp(2*(1 + k - nu),1) + (-1 + n - q)*psiexp(2*(1 + Max(0,-j + k)),1) + Max(0,-j + k)*(-1 + psiexp(2,1))*psiexp(2*(1 + Max(0,-j + k)),1) + (-n + q)*psiexp(2*(2 + Max(0,-j + k)),1)) - vEG2G2G2G2cAWBDCH*b[kappa]*Power2(1/(psiexp(1,2) - psiexp(1,1)*psiexp(1,3)))*psiexp(j - k - nu,3)*psiexp(-kappa + nu,1)*psiexp(-nu - Max(0,-j + k),2)*(-(Max(0,-j + k)*(-psiexp(1,2) + psiexp(1,1)*psiexp(1,3))*psiexp(k,2)*psiexp(1 + Max(0,-j + k),1)*psiexp(1 + nu + Max(0,-j + k),3)) + psiexp(-nu,1)*((-1 + k - n - nu + q)*psiexp(1 + k,1)*psiexp(1 + k,3)*psiexp(1 + nu + Max(0,-j + k),2) + (1 - n + q)*psiexp(1 + k,2)*psiexp(1 + nu + Max(0,-j + k),1)*psiexp(1 + nu + Max(0,-j + k),3) + (-k + n + nu - q)*psiexp(k,1)*psiexp(k,3)*psiexp(2 + nu + Max(0,-j + k),2) + (n - q)*psiexp(k,2)*psiexp(2 + nu + Max(0,-j + k),1)*psiexp(2 + nu + Max(0,-j + k),3))) + vEG2G2G2G2cAWAHBD*b[kappa]*Power2(1/(-psiexp(1,2) + psiexp(2,1)))*psiexp(j - k - kappa + nu,1)*psiexp(-Max(0,-j + k),2)*((1 - k + n + nu - q)*psiexp(2*(1 + k - nu),1)*psiexp(1 + Max(0,-j + k),2) + Max(0,-j + k)*(-psiexp(1,2) + psiexp(2,1))*psiexp(k - nu,2)*psiexp(2*(1 + Max(0,-j + k)),1) + (-1 + n - q)*psiexp(1 + k - nu,2)*psiexp(2*(1 + Max(0,-j + k)),1) + (k - n - nu + q)*psiexp(2*(k - nu),1)*psiexp(2 + Max(0,-j + k),2) + (-n + q)*psiexp(k - nu,2)*psiexp(4 + 2*Max(0,-j + k),1));
					}
					sumInuP1 = sumIkappa;
				}
				sumInu += sumInuP1*b[nu];
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += MMjkP1/(n - q);
	}
}

for (int j = 3; j <= q; j++)
{
	for (int k = 2; k <= q; k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = 0; nu <= Min(-3 + j,-2 + k); nu++)
			{
				scalar sumInuP1;
				{
					scalar sumIkappa = 0;
					for (int kappa = Max(0,1 - j + k + nu); kappa <= -2 + k; kappa++)
					{
						sumIkappa += -(vEG2G2G2G2c*(-1 + k - nu - Max(0,kappa - nu))*(k - 2*n - nu + 2*q + Max(0,kappa - nu))*b[kappa])/2. - (vEG2G2G2G2cAWAD*(-1 + k - nu - Max(0,kappa - nu))*(k - 2*n - nu + 2*q + Max(0,kappa - nu))*b[kappa]*psiexp(k - kappa,1))/2. - (vEG2G2G2G2cADAH*(-k + 2*n + nu - 2*q - Max(0,kappa - nu))*(1 - k + nu + Max(0,kappa - nu))*b[kappa]*psiexp(j - nu,1))/2. + vEG2G2G2G2cAH*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(j - k,1)*((k - n - nu + q)*psiexp(k - nu,1) + (1 - k + n + nu - q)*psiexp(1 + k - nu,1) + (-1 + n - q)*psiexp(1 + Max(0,kappa - nu),1) + Max(0,kappa - nu)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,kappa - nu),1) + (-n + q)*psiexp(2 + Max(0,kappa - nu),1)) + vEG2G2G2G2cAWADAH*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(j - kappa,1)*((k - n - nu + q)*psiexp(k - nu,1) + (1 - k + n + nu - q)*psiexp(1 + k - nu,1) + (-1 + n - q)*psiexp(1 + Max(0,kappa - nu),1) + Max(0,kappa - nu)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,kappa - nu),1) + (-n + q)*psiexp(2 + Max(0,kappa - nu),1)) + vEG2G2G2G2cAW*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-kappa + nu,1)*((k - n - nu + q)*psiexp(k - nu,1) + (1 - k + n + nu - q)*psiexp(1 + k - nu,1) + (-1 + n - q)*psiexp(1 + Max(0,kappa - nu),1) + Max(0,kappa - nu)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,kappa - nu),1) + (-n + q)*psiexp(2 + Max(0,kappa - nu),1)) + vEG2G2G2G2cAWBDBH*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(j - nu,2)*psiexp(-kappa + nu,1)*((k - n - nu + q)*psiexp(k - nu,1) + (1 - k + n + nu - q)*psiexp(1 + k - nu,1) + (-1 + n - q)*psiexp(1 + Max(0,kappa - nu),1) + Max(0,kappa - nu)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,kappa - nu),1) + (-n + q)*psiexp(2 + Max(0,kappa - nu),1)) + vEG2G2G2G2cAD*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-Max(0,kappa - nu),1)*((-n + q)*psiexp(k - nu,1) - Max(0,kappa - nu)*(-1 + psiexp(1,1))*psiexp(k - nu,1) + (-1 + n - q)*psiexp(1 + k - nu,1) + (1 - k + n + nu - q)*psiexp(1 + Max(0,kappa - nu),1) + (k - n - nu + q)*psiexp(2 + Max(0,kappa - nu),1)) + vEG2G2G2G2cAWADBH*b[kappa]*Power2(1/(-1 + psiexp(1,2)))*psiexp(j - k,2)*psiexp(k - kappa,1)*((k - n - nu + q)*psiexp(k - nu,2) + (1 - k + n + nu - q)*psiexp(1 + k - nu,2) + (-1 + n - q)*psiexp(1 + Max(0,kappa - nu),2) + Max(0,kappa - nu)*(-1 + psiexp(1,2))*psiexp(1 + Max(0,kappa - nu),2) + (-n + q)*psiexp(2 + Max(0,kappa - nu),2)) - vEG2G2G2G2cAWBD*b[kappa]*Power2(1/(psiexp(1,1) - psiexp(1,2)))*psiexp(-kappa + nu,1)*psiexp(-Max(0,kappa - nu),2)*(-(Max(0,kappa - nu)*(psiexp(1,1) - psiexp(1,2))*psiexp(k - nu,2)*psiexp(1 + Max(0,kappa - nu),1)) + (1 - n + q)*psiexp(1 + k - nu,2)*psiexp(1 + Max(0,kappa - nu),1) + (-1 + k - n - nu + q)*psiexp(1 + k - nu,1)*psiexp(1 + Max(0,kappa - nu),2) + (n - q)*psiexp(k - nu,2)*psiexp(2 + Max(0,kappa - nu),1) + (-k + n + nu - q)*psiexp(k - nu,1)*psiexp(2 + Max(0,kappa - nu),2)) + vEG2G2G2G2cADBH*b[kappa]*Power2(1/(psiexp(1,1) - psiexp(1,2)))*psiexp(j - k,2)*psiexp(-Max(0,kappa - nu),1)*((1 - k + n + nu - q)*psiexp(1 + k - nu,2)*psiexp(1 + Max(0,kappa - nu),1) + Max(0,kappa - nu)*(-psiexp(1,1) + psiexp(1,2))*psiexp(k - nu,1)*psiexp(1 + Max(0,kappa - nu),2) + (-1 + n - q)*psiexp(1 + k - nu,1)*psiexp(1 + Max(0,kappa - nu),2) + (k - n - nu + q)*psiexp(k - nu,2)*psiexp(2 + Max(0,kappa - nu),1) + (-n + q)*psiexp(k - nu,1)*psiexp(2 + Max(0,kappa - nu),2)) + vEG2G2G2G2cAWAH*b[kappa]*Power2(1/(-1 + psiexp(2,1)))*psiexp(j - k - kappa + nu,1)*((k - n - nu + q)*psiexp(2*(k - nu),1) + (1 - k + n + nu - q)*psiexp(2*(1 + k - nu),1) + (-1 + n - q)*psiexp(2*(1 + Max(0,kappa - nu)),1) + Max(0,kappa - nu)*(-1 + psiexp(2,1))*psiexp(2*(1 + Max(0,kappa - nu)),1) + (-n + q)*psiexp(2*(2 + Max(0,kappa - nu)),1)) - vEG2G2G2G2cAWBDCH*b[kappa]*Power2(1/(psiexp(1,2) - psiexp(1,1)*psiexp(1,3)))*psiexp(j - k - nu,3)*psiexp(-kappa + nu,1)*psiexp(-nu - Max(0,kappa - nu),2)*(-(Max(0,kappa - nu)*(-psiexp(1,2) + psiexp(1,1)*psiexp(1,3))*psiexp(k,2)*psiexp(1 + Max(0,kappa - nu),1)*psiexp(1 + nu + Max(0,kappa - nu),3)) + psiexp(-nu,1)*((-1 + k - n - nu + q)*psiexp(1 + k,1)*psiexp(1 + k,3)*psiexp(1 + nu + Max(0,kappa - nu),2) + (1 - n + q)*psiexp(1 + k,2)*psiexp(1 + nu + Max(0,kappa - nu),1)*psiexp(1 + nu + Max(0,kappa - nu),3) + (-k + n + nu - q)*psiexp(k,1)*psiexp(k,3)*psiexp(2 + nu + Max(0,kappa - nu),2) + (n - q)*psiexp(k,2)*psiexp(2 + nu + Max(0,kappa - nu),1)*psiexp(2 + nu + Max(0,kappa - nu),3))) + vEG2G2G2G2cAWAHBD*b[kappa]*Power2(1/(-psiexp(1,2) + psiexp(2,1)))*psiexp(j - k - kappa + nu,1)*psiexp(-Max(0,kappa - nu),2)*((1 - k + n + nu - q)*psiexp(2*(1 + k - nu),1)*psiexp(1 + Max(0,kappa - nu),2) + Max(0,kappa - nu)*(-psiexp(1,2) + psiexp(2,1))*psiexp(k - nu,2)*psiexp(2*(1 + Max(0,kappa - nu)),1) + (-1 + n - q)*psiexp(1 + k - nu,2)*psiexp(2*(1 + Max(0,kappa - nu)),1) + (k - n - nu + q)*psiexp(2*(k - nu),1)*psiexp(2 + Max(0,kappa - nu),2) + (-n + q)*psiexp(k - nu,2)*psiexp(4 + 2*Max(0,kappa - nu),1));
					}
					sumInuP1 = sumIkappa;
				}
				sumInu += sumInuP1*b[nu];
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += MMjkP1/(n - q);
	}
}

for (int j = 2; j <= -1 + q; j++)
{
	for (int k = 2; k <= -1 + q; k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = 0; nu <= Min(-2 + j,-2 + k); nu++)
			{
				scalar sumInuP1;
				{
					scalar sumIkappa = 0;
					for (int kappa = Max(1 + j,1 + k); kappa <= Min(-1 + j + k - nu,q); kappa++)
					{
						sumIkappa += (j - kappa + n - q)*vEG2G2G4c*b[kappa] + (j - kappa + n - q)*vEG2G2G4cAD*b[kappa]*psiexp(-k + kappa,1) + (j - kappa + n - q)*vEG2G2G4cAWAD*b[kappa]*psiexp(j - nu,1) + (j - kappa + n - q)*vEG2G2G4cAW*b[kappa]*psiexp(j + k - kappa - nu,1) + (j - kappa + n - q)*vEG2G2G4cAWBD*b[kappa]*psiexp(-k + kappa,2)*psiexp(j + k - kappa - nu,1);
					}
					sumInuP1 = sumIkappa;
				}
				sumInu += sumInuP1*b[nu];
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += MMjkP1/(n - q);
	}
}

for (int j = 2; j <= q; j++)
{
	for (int k = 2; k <= q; k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = 0; nu <= Min(-2 + j,-2 + k); nu++)
			{
				sumInu += -(vEG2G4G2c*(-1 + k - nu - Max(0,-j + k))*(k - 2*n - nu + 2*q + Max(0,-j + k))*b[nu])/2. - (vEG2G4G2cAWAH*(-k + 2*n + nu - 2*q - Max(0,-j + k))*(1 - k + nu + Max(0,-j + k))*b[nu]*psiexp(j - nu,1))/2. + vEG2G4G2cAH*b[nu]*Power2(1/(-1 + psiexp(1,1)))*psiexp(j - k,1)*((k - n - nu + q)*psiexp(k - nu,1) + (1 - k + n + nu - q)*psiexp(1 + k - nu,1) + (-1 + n - q)*psiexp(1 + Max(0,-j + k),1) + Max(0,-j + k)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,-j + k),1) + (-n + q)*psiexp(2 + Max(0,-j + k),1)) + vEG2G4G2cAW*b[nu]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-Max(0,-j + k),1)*((-n + q)*psiexp(k - nu,1) - Max(0,-j + k)*(-1 + psiexp(1,1))*psiexp(k - nu,1) + (-1 + n - q)*psiexp(1 + k - nu,1) + (1 - k + n + nu - q)*psiexp(1 + Max(0,-j + k),1) + (k - n - nu + q)*psiexp(2 + Max(0,-j + k),1)) + vEG2G4G2cBH*b[nu]*Power2(1/(-1 + psiexp(1,2)))*psiexp(j - k,2)*((k - n - nu + q)*psiexp(k - nu,2) + (1 - k + n + nu - q)*psiexp(1 + k - nu,2) + (-1 + n - q)*psiexp(1 + Max(0,-j + k),2) + Max(0,-j + k)*(-1 + psiexp(1,2))*psiexp(1 + Max(0,-j + k),2) + (-n + q)*psiexp(2 + Max(0,-j + k),2)) + vEG2G4G2cAWBH*b[nu]*Power2(1/(psiexp(1,1) - psiexp(1,2)))*psiexp(j - k,2)*psiexp(-Max(0,-j + k),1)*((1 - k + n + nu - q)*psiexp(1 + k - nu,2)*psiexp(1 + Max(0,-j + k),1) + Max(0,-j + k)*(-psiexp(1,1) + psiexp(1,2))*psiexp(k - nu,1)*psiexp(1 + Max(0,-j + k),2) + (-1 + n - q)*psiexp(1 + k - nu,1)*psiexp(1 + Max(0,-j + k),2) + (k - n - nu + q)*psiexp(k - nu,2)*psiexp(2 + Max(0,-j + k),1) + (-n + q)*psiexp(k - nu,1)*psiexp(2 + Max(0,-j + k),2)) + vEG2G4G2cAWCH*b[nu]*Power2(1/(psiexp(1,1) - psiexp(1,3)))*psiexp(j - k,3)*psiexp(-Max(0,-j + k),1)*((1 - k + n + nu - q)*psiexp(1 + k - nu,3)*psiexp(1 + Max(0,-j + k),1) + Max(0,-j + k)*(-psiexp(1,1) + psiexp(1,3))*psiexp(k - nu,1)*psiexp(1 + Max(0,-j + k),3) + (-1 + n - q)*psiexp(1 + k - nu,1)*psiexp(1 + Max(0,-j + k),3) + (k - n - nu + q)*psiexp(k - nu,3)*psiexp(2 + Max(0,-j + k),1) + (-n + q)*psiexp(k - nu,1)*psiexp(2 + Max(0,-j + k),3));
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += (MMjkP1*b[k])/(n - q);
	}
}

for (int j = 2; j <= q; j++)
{
	for (int k = 2; k <= Min(-2 + j + n - q,q); k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = 0; nu <= Min(-2 + j,-2 + k); nu++)
			{
				scalar sumInuP1;
				{
					scalar sumIkappa = 0;
					for (int kappa = Max3(0,1 + nu,1 - j + k + nu); kappa <= Min(-1 + k,-1 + n + nu - q); kappa++)
					{
						sumIkappa += (-kappa + n + nu - q)*vEG4G2G2c*b[kappa] - (kappa - n - nu + q)*vEG4G2G2cAD*b[kappa]*psiexp(k - kappa,1) - (kappa - n - nu + q)*vEG4G2G2cBD*b[kappa]*psiexp(k - kappa,2) + (-kappa + n + nu - q)*vEG4G2G2cADAH*b[kappa]*psiexp(j - nu,1) + (-kappa + n + nu - q)*vEG4G2G2cBDBH*b[kappa]*psiexp(j - nu,2) - (kappa - n - nu + q)*vEG4G2G2cAH*b[kappa]*psiexp(j - k + kappa - nu,1) + (-kappa + n + nu - q)*vEG4G2G2cAHBD*b[kappa]*psiexp(k - kappa,2)*psiexp(j - k + kappa - nu,1) + (-kappa + n + nu - q)*vEG4G2G2cADBH*b[kappa]*psiexp(k - kappa,1)*psiexp(j - k + kappa - nu,2) + (-kappa + n + nu - q)*vEG4G2G2cBDCH*b[kappa]*psiexp(k - kappa,2)*psiexp(j - k + kappa - nu,3);
					}
					sumInuP1 = sumIkappa;
				}
				sumInu += sumInuP1*b[nu];
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += MMjkP1/(n - q);
	}
}

for (int j = 1; j <= -2 + q; j++)
{
	for (int k = 1; k <= -3 + q; k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = Max(0,2 + j + k - q); nu <= -1 + j; nu++)
			{
				scalar sumInuP1;
				{
					scalar sumIkappa = 0;
					for (int kappa = Max(2 + j,2 + j + k - nu); kappa <= q; kappa++)
					{
						sumIkappa += -(vEG2G2G2G2c*(j - kappa + 2*n - 2*q - Max(0,k - nu))*(1 + j - kappa + Max(0,k - nu))*b[kappa])/2. - (vEG2G2G2G2cAWADAH*(j - kappa + 2*n - 2*q - Max(0,k - nu))*(1 + j - kappa + Max(0,k - nu))*b[kappa]*psiexp(-k + kappa,1))/2. - (vEG2G2G2G2cAD*(j - kappa + 2*n - 2*q - Max(0,k - nu))*(1 + j - kappa + Max(0,k - nu))*b[kappa]*psiexp(j - nu,1))/2. - (vEG2G2G2G2cAWAH*(j - kappa + 2*n - 2*q - Max(0,k - nu))*(1 + j - kappa + Max(0,k - nu))*b[kappa]*psiexp(-j - k + kappa + nu,1))/2. - (vEG2G2G2G2cAWAHBD*(j - kappa + 2*n - 2*q - Max(0,k - nu))*(1 + j - kappa + Max(0,k - nu))*b[kappa]*psiexp(j - nu,2)*psiexp(-j - k + kappa + nu,1))/2. - vEG2G2G2G2cAH*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-Max(0,k - nu),1)*((n - q)*psiexp(-j + kappa,1) + Max(0,k - nu)*(-1 + psiexp(1,1))*psiexp(-j + kappa,1) + (1 - n + q)*psiexp(1 - j + kappa,1) + (-1 - j + kappa - n + q)*psiexp(1 + Max(0,k - nu),1) + (j - kappa + n - q)*psiexp(2 + Max(0,k - nu),1)) + vEG2G2G2G2cAWAD*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(j - k,1)*((-j + kappa - n + q)*psiexp(-j + kappa,1) + (1 + j - kappa + n - q)*psiexp(1 - j + kappa,1) + (-1 + n - q)*psiexp(1 + Max(0,k - nu),1) + Max(0,k - nu)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,k - nu),1) + (-n + q)*psiexp(2 + Max(0,k - nu),1)) + vEG2G2G2G2cAW*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-k + nu,1)*((-j + kappa - n + q)*psiexp(-j + kappa,1) + (1 + j - kappa + n - q)*psiexp(1 - j + kappa,1) + (-1 + n - q)*psiexp(1 + Max(0,k - nu),1) + Max(0,k - nu)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,k - nu),1) + (-n + q)*psiexp(2 + Max(0,k - nu),1)) + vEG2G2G2G2cAWBD*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(j - nu,2)*psiexp(-k + nu,1)*((-j + kappa - n + q)*psiexp(-j + kappa,1) + (1 + j - kappa + n - q)*psiexp(1 - j + kappa,1) + (-1 + n - q)*psiexp(1 + Max(0,k - nu),1) + Max(0,k - nu)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,k - nu),1) + (-n + q)*psiexp(2 + Max(0,k - nu),1)) + vEG2G2G2G2cADAH*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(j - nu - Max(0,k - nu),1)*((-n + q)*psiexp(-j + kappa,1) - Max(0,k - nu)*(-1 + psiexp(1,1))*psiexp(-j + kappa,1) + (-1 + n - q)*psiexp(1 - j + kappa,1) + (1 + j - kappa + n - q)*psiexp(1 + Max(0,k - nu),1) + (-j + kappa - n + q)*psiexp(2 + Max(0,k - nu),1)) + vEG2G2G2G2cADBH*b[kappa]*Power2(1/(-1 + psiexp(1,2)))*psiexp(j - nu,1)*psiexp(-Max(0,k - nu),2)*((-n + q)*psiexp(-j + kappa,2) - Max(0,k - nu)*(-1 + psiexp(1,2))*psiexp(-j + kappa,2) + (-1 + n - q)*psiexp(1 - j + kappa,2) + (1 + j - kappa + n - q)*psiexp(1 + Max(0,k - nu),2) + (-j + kappa - n + q)*psiexp(2 + Max(0,k - nu),2)) - vEG2G2G2G2cAWADBH*b[kappa]*Power2(1/(psiexp(1,1) - psiexp(1,2)))*psiexp(j - k,1)*psiexp(-Max(0,k - nu),2)*(-(Max(0,k - nu)*(psiexp(1,1) - psiexp(1,2))*psiexp(-j + kappa,2)*psiexp(1 + Max(0,k - nu),1)) + (1 - n + q)*psiexp(1 - j + kappa,2)*psiexp(1 + Max(0,k - nu),1) + (-1 - j + kappa - n + q)*psiexp(1 - j + kappa,1)*psiexp(1 + Max(0,k - nu),2) + (n - q)*psiexp(-j + kappa,2)*psiexp(2 + Max(0,k - nu),1) + (j - kappa + n - q)*psiexp(-j + kappa,1)*psiexp(2 + Max(0,k - nu),2)) + vEG2G2G2G2cAWBDBH*b[kappa]*Power2(1/(psiexp(1,1) - psiexp(1,2)))*psiexp(-k + nu,1)*psiexp(j - nu - Max(0,k - nu),2)*(Max(0,k - nu)*(psiexp(1,1) - psiexp(1,2))*psiexp(-j + kappa,2)*psiexp(1 + Max(0,k - nu),1) + (-1 + n - q)*psiexp(1 - j + kappa,2)*psiexp(1 + Max(0,k - nu),1) + (1 + j - kappa + n - q)*psiexp(1 - j + kappa,1)*psiexp(1 + Max(0,k - nu),2) + (-n + q)*psiexp(-j + kappa,2)*psiexp(2 + Max(0,k - nu),1) + (-j + kappa - n + q)*psiexp(-j + kappa,1)*psiexp(2 + Max(0,k - nu),2)) - vEG2G2G2G2cAWBDCH*b[kappa]*Power2(1/(psiexp(1,1) - psiexp(1,3)))*psiexp(j - nu,2)*psiexp(-k + nu,1)*psiexp(-Max(0,k - nu),3)*(-(Max(0,k - nu)*(psiexp(1,1) - psiexp(1,3))*psiexp(-j + kappa,3)*psiexp(1 + Max(0,k - nu),1)) + (1 - n + q)*psiexp(1 - j + kappa,3)*psiexp(1 + Max(0,k - nu),1) + (-1 - j + kappa - n + q)*psiexp(1 - j + kappa,1)*psiexp(1 + Max(0,k - nu),3) + (n - q)*psiexp(-j + kappa,3)*psiexp(2 + Max(0,k - nu),1) + (j - kappa + n - q)*psiexp(-j + kappa,1)*psiexp(2 + Max(0,k - nu),3));
					}
					sumInuP1 = sumIkappa;
				}
				sumInu += sumInuP1*b[nu];
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += MMjkP1/(n - q);
	}
}

for (int j = 2; j <= q; j++)
{
	for (int k = 1; k <= Min(-3 + q,-3 - j + 2*q); k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = Max(0,1 + j + k - q); nu <= Min(-2 + j,-2 + q); nu++)
			{
				scalar sumInuP1;
				{
					scalar sumIkappa = 0;
					for (int kappa = Max(1 + j + k - nu,2 + nu); kappa <= q; kappa++)
					{
						sumIkappa += -(vEG2G2G2G2c*(-1 + kappa - nu - Max(0,-j + kappa))*(kappa - 2*n - nu + 2*q + Max(0,-j + kappa))*b[kappa])/2. - (vEG2G2G2G2cAWAD*(-1 + kappa - nu - Max(0,-j + kappa))*(kappa - 2*n - nu + 2*q + Max(0,-j + kappa))*b[kappa]*psiexp(-k + kappa,1))/2. - (vEG2G2G2G2cADAH*(-kappa + 2*n + nu - 2*q - Max(0,-j + kappa))*(1 - kappa + nu + Max(0,-j + kappa))*b[kappa]*psiexp(j - nu,1))/2. + vEG2G2G2G2cAWADAH*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(j - k,1)*((kappa - n - nu + q)*psiexp(kappa - nu,1) + (1 - kappa + n + nu - q)*psiexp(1 + kappa - nu,1) + (-1 + n - q)*psiexp(1 + Max(0,-j + kappa),1) + Max(0,-j + kappa)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,-j + kappa),1) + (-n + q)*psiexp(2 + Max(0,-j + kappa),1)) + vEG2G2G2G2cAH*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(j - kappa,1)*((kappa - n - nu + q)*psiexp(kappa - nu,1) + (1 - kappa + n + nu - q)*psiexp(1 + kappa - nu,1) + (-1 + n - q)*psiexp(1 + Max(0,-j + kappa),1) + Max(0,-j + kappa)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,-j + kappa),1) + (-n + q)*psiexp(2 + Max(0,-j + kappa),1)) + vEG2G2G2G2cAW*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-k + nu,1)*((kappa - n - nu + q)*psiexp(kappa - nu,1) + (1 - kappa + n + nu - q)*psiexp(1 + kappa - nu,1) + (-1 + n - q)*psiexp(1 + Max(0,-j + kappa),1) + Max(0,-j + kappa)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,-j + kappa),1) + (-n + q)*psiexp(2 + Max(0,-j + kappa),1)) + vEG2G2G2G2cAWBDBH*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(j - nu,2)*psiexp(-k + nu,1)*((kappa - n - nu + q)*psiexp(kappa - nu,1) + (1 - kappa + n + nu - q)*psiexp(1 + kappa - nu,1) + (-1 + n - q)*psiexp(1 + Max(0,-j + kappa),1) + Max(0,-j + kappa)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,-j + kappa),1) + (-n + q)*psiexp(2 + Max(0,-j + kappa),1)) + vEG2G2G2G2cAD*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-Max(0,-j + kappa),1)*((-n + q)*psiexp(kappa - nu,1) - Max(0,-j + kappa)*(-1 + psiexp(1,1))*psiexp(kappa - nu,1) + (-1 + n - q)*psiexp(1 + kappa - nu,1) + (1 - kappa + n + nu - q)*psiexp(1 + Max(0,-j + kappa),1) + (kappa - n - nu + q)*psiexp(2 + Max(0,-j + kappa),1)) + vEG2G2G2G2cAWADBH*b[kappa]*Power2(1/(-1 + psiexp(1,2)))*psiexp(j - kappa,2)*psiexp(-k + kappa,1)*((kappa - n - nu + q)*psiexp(kappa - nu,2) + (1 - kappa + n + nu - q)*psiexp(1 + kappa - nu,2) + (-1 + n - q)*psiexp(1 + Max(0,-j + kappa),2) + Max(0,-j + kappa)*(-1 + psiexp(1,2))*psiexp(1 + Max(0,-j + kappa),2) + (-n + q)*psiexp(2 + Max(0,-j + kappa),2)) - vEG2G2G2G2cAWBD*b[kappa]*Power2(1/(psiexp(1,1) - psiexp(1,2)))*psiexp(-k + nu,1)*psiexp(-Max(0,-j + kappa),2)*(-(Max(0,-j + kappa)*(psiexp(1,1) - psiexp(1,2))*psiexp(kappa - nu,2)*psiexp(1 + Max(0,-j + kappa),1)) + (1 - n + q)*psiexp(1 + kappa - nu,2)*psiexp(1 + Max(0,-j + kappa),1) + (-1 + kappa - n - nu + q)*psiexp(1 + kappa - nu,1)*psiexp(1 + Max(0,-j + kappa),2) + (n - q)*psiexp(kappa - nu,2)*psiexp(2 + Max(0,-j + kappa),1) + (-kappa + n + nu - q)*psiexp(kappa - nu,1)*psiexp(2 + Max(0,-j + kappa),2)) + vEG2G2G2G2cADBH*b[kappa]*Power2(1/(psiexp(1,1) - psiexp(1,2)))*psiexp(j - kappa,2)*psiexp(-Max(0,-j + kappa),1)*((1 - kappa + n + nu - q)*psiexp(1 + kappa - nu,2)*psiexp(1 + Max(0,-j + kappa),1) + Max(0,-j + kappa)*(-psiexp(1,1) + psiexp(1,2))*psiexp(kappa - nu,1)*psiexp(1 + Max(0,-j + kappa),2) + (-1 + n - q)*psiexp(1 + kappa - nu,1)*psiexp(1 + Max(0,-j + kappa),2) + (kappa - n - nu + q)*psiexp(kappa - nu,2)*psiexp(2 + Max(0,-j + kappa),1) + (-n + q)*psiexp(kappa - nu,1)*psiexp(2 + Max(0,-j + kappa),2)) + vEG2G2G2G2cAWAH*b[kappa]*Power2(1/(-1 + psiexp(2,1)))*psiexp(j - k - kappa + nu,1)*((kappa - n - nu + q)*psiexp(2*(kappa - nu),1) + (1 - kappa + n + nu - q)*psiexp(2*(1 + kappa - nu),1) + (-1 + n - q)*psiexp(2*(1 + Max(0,-j + kappa)),1) + Max(0,-j + kappa)*(-1 + psiexp(2,1))*psiexp(2*(1 + Max(0,-j + kappa)),1) + (-n + q)*psiexp(2*(2 + Max(0,-j + kappa)),1)) - vEG2G2G2G2cAWBDCH*b[kappa]*Power2(1/(psiexp(1,2) - psiexp(1,1)*psiexp(1,3)))*psiexp(j - kappa - nu,3)*psiexp(-k + nu,1)*psiexp(-nu - Max(0,-j + kappa),2)*(-(Max(0,-j + kappa)*(-psiexp(1,2) + psiexp(1,1)*psiexp(1,3))*psiexp(kappa,2)*psiexp(1 + Max(0,-j + kappa),1)*psiexp(1 + nu + Max(0,-j + kappa),3)) + psiexp(-nu,1)*((-1 + kappa - n - nu + q)*psiexp(1 + kappa,1)*psiexp(1 + kappa,3)*psiexp(1 + nu + Max(0,-j + kappa),2) + (1 - n + q)*psiexp(1 + kappa,2)*psiexp(1 + nu + Max(0,-j + kappa),1)*psiexp(1 + nu + Max(0,-j + kappa),3) + (-kappa + n + nu - q)*psiexp(kappa,1)*psiexp(kappa,3)*psiexp(2 + nu + Max(0,-j + kappa),2) + (n - q)*psiexp(kappa,2)*psiexp(2 + nu + Max(0,-j + kappa),1)*psiexp(2 + nu + Max(0,-j + kappa),3))) + vEG2G2G2G2cAWAHBD*b[kappa]*Power2(1/(-psiexp(1,2) + psiexp(2,1)))*psiexp(j - k - kappa + nu,1)*psiexp(-Max(0,-j + kappa),2)*((1 - kappa + n + nu - q)*psiexp(2*(1 + kappa - nu),1)*psiexp(1 + Max(0,-j + kappa),2) + Max(0,-j + kappa)*(-psiexp(1,2) + psiexp(2,1))*psiexp(kappa - nu,2)*psiexp(2*(1 + Max(0,-j + kappa)),1) + (-1 + n - q)*psiexp(1 + kappa - nu,2)*psiexp(2*(1 + Max(0,-j + kappa)),1) + (kappa - n - nu + q)*psiexp(2*(kappa - nu),1)*psiexp(2 + Max(0,-j + kappa),2) + (-n + q)*psiexp(kappa - nu,2)*psiexp(4 + 2*Max(0,-j + kappa),1));
					}
					sumInuP1 = sumIkappa;
				}
				sumInu += sumInuP1*b[nu];
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += MMjkP1/(n - q);
	}
}

for (int j = 2; j <= q; j++)
{
	for (int k = Max(1,2 - j); k <= -2 + q; k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = 0; nu <= Floor(Min3(-2 + j,(-2 + j + k)/2.,-2 + q)); nu++)
			{
				scalar sumInuP1;
				{
					scalar sumIkappa = 0;
					for (int kappa = Max(2 + k,2 + nu); kappa <= Min(j + k - nu,q); kappa++)
					{
						sumIkappa += -(vEG2G2G2G2c*(-1 + kappa - nu - Max(0,k - nu))*(kappa - 2*n - nu + 2*q + Max(0,k - nu))*b[kappa])/2. - (vEG2G2G2G2cAWAD*(-1 + kappa - nu - Max(0,k - nu))*(kappa - 2*n - nu + 2*q + Max(0,k - nu))*b[kappa]*psiexp(-k + kappa,1))/2. - (vEG2G2G2G2cADAH*(-kappa + 2*n + nu - 2*q - Max(0,k - nu))*(1 - kappa + nu + Max(0,k - nu))*b[kappa]*psiexp(j - nu,1))/2. + vEG2G2G2G2cAWADAH*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(j - k,1)*((kappa - n - nu + q)*psiexp(kappa - nu,1) + (1 - kappa + n + nu - q)*psiexp(1 + kappa - nu,1) + (-1 + n - q)*psiexp(1 + Max(0,k - nu),1) + Max(0,k - nu)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,k - nu),1) + (-n + q)*psiexp(2 + Max(0,k - nu),1)) + vEG2G2G2G2cAH*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(j - kappa,1)*((kappa - n - nu + q)*psiexp(kappa - nu,1) + (1 - kappa + n + nu - q)*psiexp(1 + kappa - nu,1) + (-1 + n - q)*psiexp(1 + Max(0,k - nu),1) + Max(0,k - nu)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,k - nu),1) + (-n + q)*psiexp(2 + Max(0,k - nu),1)) + vEG2G2G2G2cAW*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-k + nu,1)*((kappa - n - nu + q)*psiexp(kappa - nu,1) + (1 - kappa + n + nu - q)*psiexp(1 + kappa - nu,1) + (-1 + n - q)*psiexp(1 + Max(0,k - nu),1) + Max(0,k - nu)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,k - nu),1) + (-n + q)*psiexp(2 + Max(0,k - nu),1)) + vEG2G2G2G2cAWBDBH*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(j - nu,2)*psiexp(-k + nu,1)*((kappa - n - nu + q)*psiexp(kappa - nu,1) + (1 - kappa + n + nu - q)*psiexp(1 + kappa - nu,1) + (-1 + n - q)*psiexp(1 + Max(0,k - nu),1) + Max(0,k - nu)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,k - nu),1) + (-n + q)*psiexp(2 + Max(0,k - nu),1)) + vEG2G2G2G2cAD*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-Max(0,k - nu),1)*((-n + q)*psiexp(kappa - nu,1) - Max(0,k - nu)*(-1 + psiexp(1,1))*psiexp(kappa - nu,1) + (-1 + n - q)*psiexp(1 + kappa - nu,1) + (1 - kappa + n + nu - q)*psiexp(1 + Max(0,k - nu),1) + (kappa - n - nu + q)*psiexp(2 + Max(0,k - nu),1)) + vEG2G2G2G2cAWADBH*b[kappa]*Power2(1/(-1 + psiexp(1,2)))*psiexp(j - kappa,2)*psiexp(-k + kappa,1)*((kappa - n - nu + q)*psiexp(kappa - nu,2) + (1 - kappa + n + nu - q)*psiexp(1 + kappa - nu,2) + (-1 + n - q)*psiexp(1 + Max(0,k - nu),2) + Max(0,k - nu)*(-1 + psiexp(1,2))*psiexp(1 + Max(0,k - nu),2) + (-n + q)*psiexp(2 + Max(0,k - nu),2)) - vEG2G2G2G2cAWBD*b[kappa]*Power2(1/(psiexp(1,1) - psiexp(1,2)))*psiexp(-k + nu,1)*psiexp(-Max(0,k - nu),2)*(-(Max(0,k - nu)*(psiexp(1,1) - psiexp(1,2))*psiexp(kappa - nu,2)*psiexp(1 + Max(0,k - nu),1)) + (1 - n + q)*psiexp(1 + kappa - nu,2)*psiexp(1 + Max(0,k - nu),1) + (-1 + kappa - n - nu + q)*psiexp(1 + kappa - nu,1)*psiexp(1 + Max(0,k - nu),2) + (n - q)*psiexp(kappa - nu,2)*psiexp(2 + Max(0,k - nu),1) + (-kappa + n + nu - q)*psiexp(kappa - nu,1)*psiexp(2 + Max(0,k - nu),2)) + vEG2G2G2G2cADBH*b[kappa]*Power2(1/(psiexp(1,1) - psiexp(1,2)))*psiexp(j - kappa,2)*psiexp(-Max(0,k - nu),1)*((1 - kappa + n + nu - q)*psiexp(1 + kappa - nu,2)*psiexp(1 + Max(0,k - nu),1) + Max(0,k - nu)*(-psiexp(1,1) + psiexp(1,2))*psiexp(kappa - nu,1)*psiexp(1 + Max(0,k - nu),2) + (-1 + n - q)*psiexp(1 + kappa - nu,1)*psiexp(1 + Max(0,k - nu),2) + (kappa - n - nu + q)*psiexp(kappa - nu,2)*psiexp(2 + Max(0,k - nu),1) + (-n + q)*psiexp(kappa - nu,1)*psiexp(2 + Max(0,k - nu),2)) + vEG2G2G2G2cAWAH*b[kappa]*Power2(1/(-1 + psiexp(2,1)))*psiexp(j - k - kappa + nu,1)*((kappa - n - nu + q)*psiexp(2*(kappa - nu),1) + (1 - kappa + n + nu - q)*psiexp(2*(1 + kappa - nu),1) + (-1 + n - q)*psiexp(2*(1 + Max(0,k - nu)),1) + Max(0,k - nu)*(-1 + psiexp(2,1))*psiexp(2*(1 + Max(0,k - nu)),1) + (-n + q)*psiexp(2*(2 + Max(0,k - nu)),1)) - vEG2G2G2G2cAWBDCH*b[kappa]*Power2(1/(psiexp(1,2) - psiexp(1,1)*psiexp(1,3)))*psiexp(j - kappa - nu,3)*psiexp(-k + nu,1)*psiexp(-nu - Max(0,k - nu),2)*(-(Max(0,k - nu)*(-psiexp(1,2) + psiexp(1,1)*psiexp(1,3))*psiexp(kappa,2)*psiexp(1 + Max(0,k - nu),1)*psiexp(1 + nu + Max(0,k - nu),3)) + psiexp(-nu,1)*((-1 + kappa - n - nu + q)*psiexp(1 + kappa,1)*psiexp(1 + kappa,3)*psiexp(1 + nu + Max(0,k - nu),2) + (1 - n + q)*psiexp(1 + kappa,2)*psiexp(1 + nu + Max(0,k - nu),1)*psiexp(1 + nu + Max(0,k - nu),3) + (-kappa + n + nu - q)*psiexp(kappa,1)*psiexp(kappa,3)*psiexp(2 + nu + Max(0,k - nu),2) + (n - q)*psiexp(kappa,2)*psiexp(2 + nu + Max(0,k - nu),1)*psiexp(2 + nu + Max(0,k - nu),3))) + vEG2G2G2G2cAWAHBD*b[kappa]*Power2(1/(-psiexp(1,2) + psiexp(2,1)))*psiexp(j - k - kappa + nu,1)*psiexp(-Max(0,k - nu),2)*((1 - kappa + n + nu - q)*psiexp(2*(1 + kappa - nu),1)*psiexp(1 + Max(0,k - nu),2) + Max(0,k - nu)*(-psiexp(1,2) + psiexp(2,1))*psiexp(kappa - nu,2)*psiexp(2*(1 + Max(0,k - nu)),1) + (-1 + n - q)*psiexp(1 + kappa - nu,2)*psiexp(2*(1 + Max(0,k - nu)),1) + (kappa - n - nu + q)*psiexp(2*(kappa - nu),1)*psiexp(2 + Max(0,k - nu),2) + (-n + q)*psiexp(kappa - nu,2)*psiexp(4 + 2*Max(0,k - nu),1));
					}
					sumInuP1 = sumIkappa;
				}
				sumInu += sumInuP1*b[nu];
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += MMjkP1/(n - q);
	}
}

for (int j = 1; j <= q; j++)
{
	for (int k = 1; k <= Min(-4 + j + n - q,-1 + q); k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = Max(0,3 + k - n + q); nu <= -1 + j; nu++)
			{
				scalar sumInuP1;
				{
					scalar sumIkappa = 0;
					for (int kappa = 1 + k; kappa <= Min(-2 + n + nu - q,q); kappa++)
					{
						sumIkappa += (vEG2G2G2G2c*b[kappa]*(-n + q - 2*n*q + (1 - 2*n + 2*q)*Max(0,kappa - nu) + Power2(n) + Power2(q) + Power2(Max(0,kappa - nu))))/2. + (vEG2G2G2G2cAW*b[kappa]*(-n + q - 2*n*q + (1 - 2*n + 2*q)*Max(0,kappa - nu) + Power2(n) + Power2(q) + Power2(Max(0,kappa - nu)))*psiexp(-k + kappa,1))/2. + (vEG2G2G2G2cAH*b[kappa]*(-n + q - 2*n*q + (1 - 2*n + 2*q)*Max(0,kappa - nu) + Power2(n) + Power2(q) + Power2(Max(0,kappa - nu)))*psiexp(j - nu,1))/2. + (vEG2G2G2G2cAWAH*b[kappa]*(-n + q - 2*n*q + (1 - 2*n + 2*q)*Max(0,kappa - nu) + Power2(n) + Power2(q) + Power2(Max(0,kappa - nu)))*psiexp(j - k + kappa - nu,1))/2. + vEG2G2G2G2cAWADAH*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(j - k,1)*(psiexp(1 + n - q,1) + (-1 + n - q)*psiexp(1 + Max(0,kappa - nu),1) + Max(0,kappa - nu)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,kappa - nu),1) + (-n + q)*psiexp(2 + Max(0,kappa - nu),1)) + vEG2G2G2G2cADAH*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(j - kappa,1)*(psiexp(1 + n - q,1) + (-1 + n - q)*psiexp(1 + Max(0,kappa - nu),1) + Max(0,kappa - nu)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,kappa - nu),1) + (-n + q)*psiexp(2 + Max(0,kappa - nu),1)) + vEG2G2G2G2cAWAD*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-k + nu,1)*(psiexp(1 + n - q,1) + (-1 + n - q)*psiexp(1 + Max(0,kappa - nu),1) + Max(0,kappa - nu)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,kappa - nu),1) + (-n + q)*psiexp(2 + Max(0,kappa - nu),1)) + vEG2G2G2G2cAWADBH*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(j - nu,2)*psiexp(-k + nu,1)*(psiexp(1 + n - q,1) + (-1 + n - q)*psiexp(1 + Max(0,kappa - nu),1) + Max(0,kappa - nu)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,kappa - nu),1) + (-n + q)*psiexp(2 + Max(0,kappa - nu),1)) + vEG2G2G2G2cAD*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-kappa + nu,1)*(psiexp(1 + n - q,1) + (-1 + n - q)*psiexp(1 + Max(0,kappa - nu),1) + Max(0,kappa - nu)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,kappa - nu),1) + (-n + q)*psiexp(2 + Max(0,kappa - nu),1)) + vEG2G2G2G2cADBH*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(j - nu,2)*psiexp(-kappa + nu,1)*(psiexp(1 + n - q,1) + (-1 + n - q)*psiexp(1 + Max(0,kappa - nu),1) + Max(0,kappa - nu)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,kappa - nu),1) + (-n + q)*psiexp(2 + Max(0,kappa - nu),1)) + vEG2G2G2G2cAWBDBH*b[kappa]*Power2(1/(-1 + psiexp(1,2)))*psiexp(j - kappa,2)*psiexp(-k + kappa,1)*(psiexp(1 + n - q,2) + (-1 + n - q)*psiexp(1 + Max(0,kappa - nu),2) + Max(0,kappa - nu)*(-1 + psiexp(1,2))*psiexp(1 + Max(0,kappa - nu),2) + (-n + q)*psiexp(2 + Max(0,kappa - nu),2)) + vEG2G2G2G2cAWBD*b[kappa]*Power2(1/(-1 + psiexp(1,2)))*psiexp(-k + kappa,1)*psiexp(-kappa + nu,2)*(psiexp(1 + n - q,2) + (-1 + n - q)*psiexp(1 + Max(0,kappa - nu),2) + Max(0,kappa - nu)*(-1 + psiexp(1,2))*psiexp(1 + Max(0,kappa - nu),2) + (-n + q)*psiexp(2 + Max(0,kappa - nu),2)) + vEG2G2G2G2cAWBDCH*b[kappa]*Power2(1/(-1 + psiexp(1,2)))*psiexp(-k + kappa,1)*psiexp(j - nu,3)*psiexp(-kappa + nu,2)*(psiexp(1 + n - q,2) + (-1 + n - q)*psiexp(1 + Max(0,kappa - nu),2) + Max(0,kappa - nu)*(-1 + psiexp(1,2))*psiexp(1 + Max(0,kappa - nu),2) + (-n + q)*psiexp(2 + Max(0,kappa - nu),2)) + vEG2G2G2G2cAWAHBD*b[kappa]*Power2(1/(-1 + psiexp(1,2)))*psiexp(j - k + kappa - nu,1)*psiexp(-kappa + nu,2)*(psiexp(1 + n - q,2) + (-1 + n - q)*psiexp(1 + Max(0,kappa - nu),2) + Max(0,kappa - nu)*(-1 + psiexp(1,2))*psiexp(1 + Max(0,kappa - nu),2) + (-n + q)*psiexp(2 + Max(0,kappa - nu),2));
					}
					sumInuP1 = sumIkappa;
				}
				sumInu += sumInuP1*b[nu];
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += MMjkP1/(n - q);
	}
}

for (int j = 1; j <= q; j++)
{
	for (int k = 1; k <= Min(-3 + j + n - q,q); k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = Max(0,2 + k - n + q); nu <= -1 + j; nu++)
			{
				scalar sumInuP1;
				{
					scalar sumIkappa = 0;
					for (int kappa = 0; kappa <= -1 + k; kappa++)
					{
						sumIkappa += (vEG2G2G2G2c*b[kappa]*(-n + q - 2*n*q + (1 - 2*n + 2*q)*Max(0,k - nu) + Power2(n) + Power2(q) + Power2(Max(0,k - nu))))/2. + (vEG2G2G2G2cAW*b[kappa]*(-n + q - 2*n*q + (1 - 2*n + 2*q)*Max(0,k - nu) + Power2(n) + Power2(q) + Power2(Max(0,k - nu)))*psiexp(k - kappa,1))/2. + (vEG2G2G2G2cAH*b[kappa]*(-n + q - 2*n*q + (1 - 2*n + 2*q)*Max(0,k - nu) + Power2(n) + Power2(q) + Power2(Max(0,k - nu)))*psiexp(j - nu,1))/2. + (vEG2G2G2G2cAWAH*b[kappa]*(-n + q - 2*n*q + (1 - 2*n + 2*q)*Max(0,k - nu) + Power2(n) + Power2(q) + Power2(Max(0,k - nu)))*psiexp(j + k - kappa - nu,1))/2. + vEG2G2G2G2cADAH*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(j - k,1)*(psiexp(1 + n - q,1) + (-1 + n - q)*psiexp(1 + Max(0,k - nu),1) + Max(0,k - nu)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,k - nu),1) + (-n + q)*psiexp(2 + Max(0,k - nu),1)) + vEG2G2G2G2cAWADAH*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(j - kappa,1)*(psiexp(1 + n - q,1) + (-1 + n - q)*psiexp(1 + Max(0,k - nu),1) + Max(0,k - nu)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,k - nu),1) + (-n + q)*psiexp(2 + Max(0,k - nu),1)) + vEG2G2G2G2cAD*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-k + nu,1)*(psiexp(1 + n - q,1) + (-1 + n - q)*psiexp(1 + Max(0,k - nu),1) + Max(0,k - nu)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,k - nu),1) + (-n + q)*psiexp(2 + Max(0,k - nu),1)) + vEG2G2G2G2cADBH*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(j - nu,2)*psiexp(-k + nu,1)*(psiexp(1 + n - q,1) + (-1 + n - q)*psiexp(1 + Max(0,k - nu),1) + Max(0,k - nu)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,k - nu),1) + (-n + q)*psiexp(2 + Max(0,k - nu),1)) + vEG2G2G2G2cAWAD*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-kappa + nu,1)*(psiexp(1 + n - q,1) + (-1 + n - q)*psiexp(1 + Max(0,k - nu),1) + Max(0,k - nu)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,k - nu),1) + (-n + q)*psiexp(2 + Max(0,k - nu),1)) + vEG2G2G2G2cAWADBH*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(j - nu,2)*psiexp(-kappa + nu,1)*(psiexp(1 + n - q,1) + (-1 + n - q)*psiexp(1 + Max(0,k - nu),1) + Max(0,k - nu)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,k - nu),1) + (-n + q)*psiexp(2 + Max(0,k - nu),1)) + vEG2G2G2G2cAWBDBH*b[kappa]*Power2(1/(-1 + psiexp(1,2)))*psiexp(j - k,2)*psiexp(k - kappa,1)*(psiexp(1 + n - q,2) + (-1 + n - q)*psiexp(1 + Max(0,k - nu),2) + Max(0,k - nu)*(-1 + psiexp(1,2))*psiexp(1 + Max(0,k - nu),2) + (-n + q)*psiexp(2 + Max(0,k - nu),2)) + vEG2G2G2G2cAWBD*b[kappa]*Power2(1/(-1 + psiexp(1,2)))*psiexp(k - kappa,1)*psiexp(-k + nu,2)*(psiexp(1 + n - q,2) + (-1 + n - q)*psiexp(1 + Max(0,k - nu),2) + Max(0,k - nu)*(-1 + psiexp(1,2))*psiexp(1 + Max(0,k - nu),2) + (-n + q)*psiexp(2 + Max(0,k - nu),2)) + vEG2G2G2G2cAWBDCH*b[kappa]*Power2(1/(-1 + psiexp(1,2)))*psiexp(k - kappa,1)*psiexp(j - nu,3)*psiexp(-k + nu,2)*(psiexp(1 + n - q,2) + (-1 + n - q)*psiexp(1 + Max(0,k - nu),2) + Max(0,k - nu)*(-1 + psiexp(1,2))*psiexp(1 + Max(0,k - nu),2) + (-n + q)*psiexp(2 + Max(0,k - nu),2)) + vEG2G2G2G2cAWAHBD*b[kappa]*Power2(1/(-1 + psiexp(1,2)))*psiexp(j + k - kappa - nu,1)*psiexp(-k + nu,2)*(psiexp(1 + n - q,2) + (-1 + n - q)*psiexp(1 + Max(0,k - nu),2) + Max(0,k - nu)*(-1 + psiexp(1,2))*psiexp(1 + Max(0,k - nu),2) + (-n + q)*psiexp(2 + Max(0,k - nu),2));
					}
					sumInuP1 = sumIkappa;
				}
				sumInu += sumInuP1*b[nu];
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += MMjkP1/(n - q);
	}
}

for (int j = 1; j <= -1 + q; j++)
{
	for (int k = 1; k <= -2 + q; k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = Max(0,1 + j + k - q); nu <= -1 + j; nu++)
			{
				scalar sumInuP1;
				{
					scalar sumIkappa = 0;
					for (int kappa = Max(1 + j,1 + j + k - nu); kappa <= q; kappa++)
					{
						sumIkappa += (j - kappa + n - q)*vEG2G2G4c*b[kappa] + (j - kappa + n - q)*vEG2G2G4cAWAD*b[kappa]*psiexp(-k + kappa,1) + (j - kappa + n - q)*vEG2G2G4cAD*b[kappa]*psiexp(j - nu,1) + (j - kappa + n - q)*vEG2G2G4cAW*b[kappa]*psiexp(-j - k + kappa + nu,1) + (j - kappa + n - q)*vEG2G2G4cAWBD*b[kappa]*psiexp(j - nu,2)*psiexp(-j - k + kappa + nu,1);
					}
					sumInuP1 = sumIkappa;
				}
				sumInu += sumInuP1*b[nu];
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += MMjkP1/(n - q);
	}
}

for (int j = 1; j <= q; j++)
{
	for (int k = 1; k <= Min(-3 + j + n - q,-1 + q); k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = Max(0,2 + k - n + q); nu <= Min(-1 + j,-1 + q); nu++)
			{
				scalar sumInuP1;
				{
					scalar sumIkappa = 0;
					for (int kappa = Max(1 + k,1 + nu); kappa <= Min(-1 + n + nu - q,q); kappa++)
					{
						sumIkappa += (-kappa + n + nu - q)*vEG2G4G2c*b[kappa] + (-kappa + n + nu - q)*vEG2G4G2cAW*b[kappa]*psiexp(-k + kappa,1) + (-kappa + n + nu - q)*vEG2G4G2cAH*b[kappa]*psiexp(j - nu,1) + (-kappa + n + nu - q)*vEG2G4G2cBH*b[kappa]*psiexp(j - nu,2) + (-kappa + n + nu - q)*vEG2G4G2cAWBH*b[kappa]*psiexp(-k + kappa,1)*psiexp(j - nu,2) + (-kappa + n + nu - q)*vEG2G4G2cAWCH*b[kappa]*psiexp(-k + kappa,1)*psiexp(j - nu,3) - (kappa - n - nu + q)*vEG2G4G2cAWAH*b[kappa]*psiexp(j - k + kappa - nu,1);
					}
					sumInuP1 = sumIkappa;
				}
				sumInu += sumInuP1*b[nu];
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += MMjkP1/(n - q);
	}
}

for (int j = 1; j <= q; j++)
{
	for (int k = 1; k <= Min(-3 + j + n - q,q); k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = Max(0,2 + k - n + q); nu <= -1 + j; nu++)
			{
				sumInu += (vEG4G2G2c*b[nu]*(-n + q - 2*n*q + (1 - 2*n + 2*q)*Max(0,k - nu) + Power2(n) + Power2(q) + Power2(Max(0,k - nu))))/2. + (vEG4G2G2cAH*b[nu]*(-n + q - 2*n*q + (1 - 2*n + 2*q)*Max(0,k - nu) + Power2(n) + Power2(q) + Power2(Max(0,k - nu)))*psiexp(j - nu,1))/2. + vEG4G2G2cADAH*b[nu]*Power2(1/(-1 + psiexp(1,1)))*psiexp(j - k,1)*(psiexp(1 + n - q,1) + (-1 + n - q)*psiexp(1 + Max(0,k - nu),1) + Max(0,k - nu)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,k - nu),1) + (-n + q)*psiexp(2 + Max(0,k - nu),1)) + vEG4G2G2cAD*b[nu]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-k + nu,1)*(psiexp(1 + n - q,1) + (-1 + n - q)*psiexp(1 + Max(0,k - nu),1) + Max(0,k - nu)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,k - nu),1) + (-n + q)*psiexp(2 + Max(0,k - nu),1)) + vEG4G2G2cADBH*b[nu]*Power2(1/(-1 + psiexp(1,1)))*psiexp(j - nu,2)*psiexp(-k + nu,1)*(psiexp(1 + n - q,1) + (-1 + n - q)*psiexp(1 + Max(0,k - nu),1) + Max(0,k - nu)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,k - nu),1) + (-n + q)*psiexp(2 + Max(0,k - nu),1)) + vEG4G2G2cBDBH*b[nu]*Power2(1/(-1 + psiexp(1,2)))*psiexp(j - k,2)*(psiexp(1 + n - q,2) + (-1 + n - q)*psiexp(1 + Max(0,k - nu),2) + Max(0,k - nu)*(-1 + psiexp(1,2))*psiexp(1 + Max(0,k - nu),2) + (-n + q)*psiexp(2 + Max(0,k - nu),2)) + vEG4G2G2cBD*b[nu]*Power2(1/(-1 + psiexp(1,2)))*psiexp(-k + nu,2)*(psiexp(1 + n - q,2) + (-1 + n - q)*psiexp(1 + Max(0,k - nu),2) + Max(0,k - nu)*(-1 + psiexp(1,2))*psiexp(1 + Max(0,k - nu),2) + (-n + q)*psiexp(2 + Max(0,k - nu),2)) + vEG4G2G2cAHBD*b[nu]*Power2(1/(-1 + psiexp(1,2)))*psiexp(j - nu,1)*psiexp(-k + nu,2)*(psiexp(1 + n - q,2) + (-1 + n - q)*psiexp(1 + Max(0,k - nu),2) + Max(0,k - nu)*(-1 + psiexp(1,2))*psiexp(1 + Max(0,k - nu),2) + (-n + q)*psiexp(2 + Max(0,k - nu),2)) + vEG4G2G2cBDCH*b[nu]*Power2(1/(-1 + psiexp(1,2)))*psiexp(j - nu,3)*psiexp(-k + nu,2)*(psiexp(1 + n - q,2) + (-1 + n - q)*psiexp(1 + Max(0,k - nu),2) + Max(0,k - nu)*(-1 + psiexp(1,2))*psiexp(1 + Max(0,k - nu),2) + (-n + q)*psiexp(2 + Max(0,k - nu),2));
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += (MMjkP1*b[k])/(n - q);
	}
}

for (int j = 1; j <= -2 + q; j++)
{
	for (int k = Max(1,1 + j); k <= -1 + q; k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = 0; nu <= -1 + j; nu++)
			{
				scalar sumInuP1;
				{
					scalar sumIkappa = 0;
					for (int kappa = 1 + k; kappa <= q; kappa++)
					{
						sumIkappa += (j - k + n - q)*vEG2G4G2c*b[kappa] + (j - k + n - q)*vEG2G4G2cAH*b[kappa]*psiexp(-k + kappa,1) + (j - k + n - q)*vEG2G4G2cBH*b[kappa]*psiexp(-k + kappa,2) + (j - k + n - q)*vEG2G4G2cAW*b[kappa]*psiexp(j - nu,1) + (j - k + n - q)*vEG2G4G2cAWBH*b[kappa]*psiexp(-k + kappa,2)*psiexp(j - nu,1) + (j - k + n - q)*vEG2G4G2cAWCH*b[kappa]*psiexp(-k + kappa,3)*psiexp(j - nu,1) + (j - k + n - q)*vEG2G4G2cAWAH*b[kappa]*psiexp(j - k + kappa - nu,1);
					}
					sumInuP1 = sumIkappa;
				}
				sumInu += sumInuP1*b[nu];
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += MMjkP1/(n - q);
	}
}

for (int j = 2; j <= -1 + q; j++)
{
	for (int k = Max(1,1 + j); k <= q; k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = 0; nu <= -2 + j; nu++)
			{
				scalar sumInuP1;
				{
					scalar sumIkappa = 0;
					for (int kappa = Max(0,1 - j + k + nu); kappa <= -1 + k; kappa++)
					{
						sumIkappa += (j - k + n - q)*vEG2G2G4c*b[kappa] + (j - k + n - q)*vEG2G2G4cAD*b[kappa]*psiexp(k - kappa,1) + (j - k + n - q)*vEG2G2G4cAWAD*b[kappa]*psiexp(j - nu,1) + (j - k + n - q)*vEG2G2G4cAW*b[kappa]*psiexp(j - k + kappa - nu,1) + (j - k + n - q)*vEG2G2G4cAWBD*b[kappa]*psiexp(k - kappa,2)*psiexp(j - k + kappa - nu,1);
					}
					sumInuP1 = sumIkappa;
				}
				sumInu += sumInuP1*b[nu];
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += MMjkP1/(n - q);
	}
}

for (int j = 1; j <= -1 + q; j++)
{
	for (int k = Max(2,1 + j); k <= q; k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = Max(0,1 + j - k); nu <= -1 + j; nu++)
			{
				scalar sumInuP1;
				{
					scalar sumIkappa = 0;
					for (int kappa = 0; kappa <= -1 - j + k + nu; kappa++)
					{
						sumIkappa += (j - k + n - q)*vEG2G2G4c*b[kappa] + (j - k + n - q)*vEG2G2G4cAWAD*b[kappa]*psiexp(k - kappa,1) + (j - k + n - q)*vEG2G2G4cAD*b[kappa]*psiexp(j - nu,1) + (j - k + n - q)*vEG2G2G4cAW*b[kappa]*psiexp(-j + k - kappa + nu,1) + (j - k + n - q)*vEG2G2G4cAWBD*b[kappa]*psiexp(j - nu,2)*psiexp(-j + k - kappa + nu,1);
					}
					sumInuP1 = sumIkappa;
				}
				sumInu += sumInuP1*b[nu];
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += MMjkP1/(n - q);
	}
}

for (int j = 1; j <= -1 + q; j++)
{
	for (int k = Max(1,1 + j); k <= q; k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = 0; nu <= -1 + j; nu++)
			{
				sumInu += vEG2G6c*b[nu] + vEG2G6cAW*b[nu]*psiexp(j - nu,1);
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += (MMjkP1*(j - k + n - q)*b[k])/(n - q);
	}
}

for (int j = 1; j <= -1 + q; j++)
{
	for (int k = Max(1,1 + j); k <= q; k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = Max(0,j - k); nu <= -1 + j; nu++)
			{
				scalar sumInuP1;
				{
					scalar sumIkappa = 0;
					for (int kappa = Max(0,-j + k + nu); kappa <= -j + k + nu; kappa++)
					{
						sumIkappa += b[kappa];
					}
					sumInuP1 = sumIkappa;
				}
				sumInu += sumInuP1*b[nu]*(vEG4G4c + vEG4G4cAD*psiexp(j - nu,1) + vEG4G4cBD*psiexp(j - nu,2));
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += (MMjkP1*(j - k + n - q))/(n - q);
	}
}

for (int j = 1; j <= -2 + q; j++)
{
	for (int k = 1; k <= Min(-2 + j + n - q,-2 + q); k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = Max3(0,1 + j + k - q,1 + k - n + q); nu <= Min(-1 + j,-1 + k); nu++)
			{
				scalar sumInuP1;
				{
					scalar sumIkappa = 0;
					for (int kappa = 1 + j + k - nu; kappa <= q; kappa++)
					{
						sumIkappa += (-k + n + nu - q)*vEG4G2G2c*b[kappa] - (k - n - nu + q)*vEG4G2G2cADAH*b[kappa]*psiexp(-k + kappa,1) - (k - n - nu + q)*vEG4G2G2cBDBH*b[kappa]*psiexp(-k + kappa,2) + (-k + n + nu - q)*vEG4G2G2cAD*b[kappa]*psiexp(j - nu,1) + (-k + n + nu - q)*vEG4G2G2cBD*b[kappa]*psiexp(j - nu,2) - (k - n - nu + q)*vEG4G2G2cAH*b[kappa]*psiexp(-j - k + kappa + nu,1) - (k - n - nu + q)*vEG4G2G2cAHBD*b[kappa]*psiexp(j - nu,2)*psiexp(-j - k + kappa + nu,1) + (-k + n + nu - q)*vEG4G2G2cADBH*b[kappa]*psiexp(j - nu,1)*psiexp(-j - k + kappa + nu,2) + (-k + n + nu - q)*vEG4G2G2cBDCH*b[kappa]*psiexp(j - nu,2)*psiexp(-j - k + kappa + nu,3);
					}
					sumInuP1 = sumIkappa;
				}
				sumInu += sumInuP1*b[nu];
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += MMjkP1/(n - q);
	}
}

for (int j = 2; j <= q; j++)
{
	for (int k = 1; k <= Min(-3 + j + n - q,-1 + q); k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = Max(0,1 + k - n + q); nu <= Min(-2 + j,-1 + k); nu++)
			{
				scalar sumInuP1;
				{
					scalar sumIkappa = 0;
					for (int kappa = 1 + k; kappa <= Min(-1 + j + k - nu,q); kappa++)
					{
						sumIkappa += (-k + n + nu - q)*vEG4G2G2c*b[kappa] - (k - n - nu + q)*vEG4G2G2cAD*b[kappa]*psiexp(-k + kappa,1) - (k - n - nu + q)*vEG4G2G2cBD*b[kappa]*psiexp(-k + kappa,2) + (-k + n + nu - q)*vEG4G2G2cADAH*b[kappa]*psiexp(j - nu,1) + (-k + n + nu - q)*vEG4G2G2cBDBH*b[kappa]*psiexp(j - nu,2) - (k - n - nu + q)*vEG4G2G2cAH*b[kappa]*psiexp(j + k - kappa - nu,1) + (-k + n + nu - q)*vEG4G2G2cAHBD*b[kappa]*psiexp(-k + kappa,2)*psiexp(j + k - kappa - nu,1) + (-k + n + nu - q)*vEG4G2G2cADBH*b[kappa]*psiexp(-k + kappa,1)*psiexp(j + k - kappa - nu,2) + (-k + n + nu - q)*vEG4G2G2cBDCH*b[kappa]*psiexp(-k + kappa,2)*psiexp(j + k - kappa - nu,3);
					}
					sumInuP1 = sumIkappa;
				}
				sumInu += sumInuP1*b[nu];
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += MMjkP1/(n - q);
	}
}

for (int j = 1; j <= q; j++)
{
	for (int k = 1; k <= Min(-2 + j + n - q,q); k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = Max(0,1 + k - n + q); nu <= Min(-1 + j,-1 + k); nu++)
			{
				scalar sumInuP1;
				{
					scalar sumIkappa = 0;
					for (int kappa = 0; kappa <= -1 + k; kappa++)
					{
						sumIkappa += (-k + n + nu - q)*vEG2G4G2c*b[kappa] - (k - n - nu + q)*vEG2G4G2cAW*b[kappa]*psiexp(k - kappa,1) + (-k + n + nu - q)*vEG2G4G2cAH*b[kappa]*psiexp(j - nu,1) + (-k + n + nu - q)*vEG2G4G2cBH*b[kappa]*psiexp(j - nu,2) - (k - n - nu + q)*vEG2G4G2cAWBH*b[kappa]*psiexp(k - kappa,1)*psiexp(j - nu,2) - (k - n - nu + q)*vEG2G4G2cAWCH*b[kappa]*psiexp(k - kappa,1)*psiexp(j - nu,3) - (k - n - nu + q)*vEG2G4G2cAWAH*b[kappa]*psiexp(j + k - kappa - nu,1);
					}
					sumInuP1 = sumIkappa;
				}
				sumInu += sumInuP1*b[nu];
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += MMjkP1/(n - q);
	}
}

for (int j = 1; j <= -1 + q; j++)
{
	for (int k = 1; k <= Min(-2 + j + n - q,-1 + q); k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = Max3(0,j + k - q,1 + k - n + q); nu <= Min(-1 + j,-1 + k); nu++)
			{
				scalar sumInuP1;
				{
					scalar sumIkappa = 0;
					for (int kappa = j + k - nu; kappa <= Min(j + k - nu,q); kappa++)
					{
						sumIkappa += b[kappa];
					}
					sumInuP1 = sumIkappa;
				}
				sumInu += (-k + n + nu - q)*sumInuP1*b[nu]*(vEG4G4c + vEG4G4cAD*psiexp(j - nu,1) + vEG4G4cBD*psiexp(j - nu,2));
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += MMjkP1/(n - q);
	}
}

for (int j = 1; j <= q; j++)
{
	for (int k = 1; k <= Min(-2 + j + n - q,q); k++)
	{
		scalar MMjkP1;
		{
			scalar sumInu = 0;
			for (int nu = Max(0,1 + k - n + q); nu <= Min(-1 + j,-1 + k); nu++)
			{
				sumInu += (-k + n + nu - q)*vEG6G2c*b[nu] + (-k + n + nu - q)*vEG6G2cAH*b[nu]*psiexp(j - nu,1) + (-k + n + nu - q)*vEG6G2cBH*b[nu]*psiexp(j - nu,2) + (-k + n + nu - q)*vEG6G2cCH*b[nu]*psiexp(j - nu,3);
			}
			MMjkP1 = sumInu;
		}
		MM(j, k) += (MMjkP1*b[k])/(n - q);
	}
}

for (int j = 1; j <= -3 + q; j++)
{
	for (int k = Max(1,2 + j); k <= -1 + q; k++)
	{
		scalar MMjkP1;
		{
			scalar sumIkappa = 0;
			for (int kappa = 1 + k; kappa <= q; kappa++)
			{
				sumIkappa += -((1 + j - k)*(j - k + 2*n - 2*q)*vEG4G2G2c*b[kappa])/2. + vEG4G2G2cAD*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*(-n + q + (-1 + n - q)*psiexp(1,1) + (1 + j - k + n - q)*psiexp(1 + j - k,1) + (-j + k - n + q)*psiexp(2 + j - k,1))*psiexp(-j + k,1) + vEG4G2G2cBD*b[kappa]*Power2(1/(-1 + psiexp(1,2)))*(-n + q + (-1 + n - q)*psiexp(1,2) + (1 + j - k + n - q)*psiexp(1 + j - k,2) + (-j + k - n + q)*psiexp(2 + j - k,2))*psiexp(-j + k,2) + vEG4G2G2cADAH*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*(-n + q + (-1 + n - q)*psiexp(1,1) + (1 + j - k + n - q)*psiexp(1 + j - k,1) + (-j + k - n + q)*psiexp(2 + j - k,1))*psiexp(-j + kappa,1) + vEG4G2G2cBDBH*b[kappa]*Power2(1/(-1 + psiexp(1,2)))*(-n + q + (-1 + n - q)*psiexp(1,2) + (1 + j - k + n - q)*psiexp(1 + j - k,2) + (-j + k - n + q)*psiexp(2 + j - k,2))*psiexp(-j + kappa,2) - ((1 + j - k)*(j - k + 2*n - 2*q)*vEG4G2G2cAH*b[kappa]*psiexp(-k + kappa,1))/2. + vEG4G2G2cAHBD*b[kappa]*Power2(1/(-1 + psiexp(1,2)))*(-n + q + (-1 + n - q)*psiexp(1,2) + (1 + j - k + n - q)*psiexp(1 + j - k,2) + (-j + k - n + q)*psiexp(2 + j - k,2))*psiexp(-j + k,2)*psiexp(-k + kappa,1) + vEG4G2G2cADBH*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*(-n + q + (-1 + n - q)*psiexp(1,1) + (1 + j - k + n - q)*psiexp(1 + j - k,1) + (-j + k - n + q)*psiexp(2 + j - k,1))*psiexp(-j + k,1)*psiexp(-k + kappa,2) + vEG4G2G2cBDCH*b[kappa]*Power2(1/(-1 + psiexp(1,2)))*(-n + q + (-1 + n - q)*psiexp(1,2) + (1 + j - k + n - q)*psiexp(1 + j - k,2) + (-j + k - n + q)*psiexp(2 + j - k,2))*psiexp(-j + k,2)*psiexp(-k + kappa,3);
			}
			MMjkP1 = sumIkappa;
		}
		MM(j, k) += (MMjkP1*b[j])/(n - q);
	}
}

for (int j = 1; j <= -3 + q; j++)
{
	for (int k = Max(1,3 + j); k <= q; k++)
	{
		scalar MMjkP1;
		{
			scalar sumIkappa = 0;
			for (int kappa = Max(0,2 + j); kappa <= -1 + k; kappa++)
			{
				sumIkappa += -((1 + j - kappa)*(j - kappa + 2*n - 2*q)*vEG4G2G2c*b[kappa])/2. + vEG4G2G2cADAH*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-j + k,1)*(-n + q + (-1 + n - q)*psiexp(1,1) + (1 + j - kappa + n - q)*psiexp(1 + j - kappa,1) + (-j + kappa - n + q)*psiexp(2 + j - kappa,1)) + vEG4G2G2cBDBH*b[kappa]*Power2(1/(-1 + psiexp(1,2)))*psiexp(-j + k,2)*(-n + q + (-1 + n - q)*psiexp(1,2) + (1 + j - kappa + n - q)*psiexp(1 + j - kappa,2) + (-j + kappa - n + q)*psiexp(2 + j - kappa,2)) - ((-1 - j + kappa)*(-j + kappa - 2*n + 2*q)*vEG4G2G2cAH*b[kappa]*psiexp(k - kappa,1))/2. + vEG4G2G2cAD*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*(-n + q + (-1 + n - q)*psiexp(1,1) + (1 + j - kappa + n - q)*psiexp(1 + j - kappa,1) + (-j + kappa - n + q)*psiexp(2 + j - kappa,1))*psiexp(-j + kappa,1) + vEG4G2G2cADBH*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*(-n + q + (-1 + n - q)*psiexp(1,1) + (1 + j - kappa + n - q)*psiexp(1 + j - kappa,1) + (-j + kappa - n + q)*psiexp(2 + j - kappa,1))*psiexp(k - kappa,2)*psiexp(-j + kappa,1) + vEG4G2G2cBD*b[kappa]*Power2(1/(-1 + psiexp(1,2)))*(-n + q + (-1 + n - q)*psiexp(1,2) + (1 + j - kappa + n - q)*psiexp(1 + j - kappa,2) + (-j + kappa - n + q)*psiexp(2 + j - kappa,2))*psiexp(-j + kappa,2) + vEG4G2G2cAHBD*b[kappa]*Power2(1/(-1 + psiexp(1,2)))*(-n + q + (-1 + n - q)*psiexp(1,2) + (1 + j - kappa + n - q)*psiexp(1 + j - kappa,2) + (-j + kappa - n + q)*psiexp(2 + j - kappa,2))*psiexp(k - kappa,1)*psiexp(-j + kappa,2) + vEG4G2G2cBDCH*b[kappa]*Power2(1/(-1 + psiexp(1,2)))*(-n + q + (-1 + n - q)*psiexp(1,2) + (1 + j - kappa + n - q)*psiexp(1 + j - kappa,2) + (-j + kappa - n + q)*psiexp(2 + j - kappa,2))*psiexp(k - kappa,3)*psiexp(-j + kappa,2);
			}
			MMjkP1 = sumIkappa;
		}
		MM(j, k) += (MMjkP1*b[j])/(n - q);
	}
}

for (int j = 1; j <= -2 + q; j++)
{
	for (int k = Max(2,2 + j); k <= q; k++)
	{
		scalar MMjkP1;
		{
			scalar sumIkappa = 0;
			for (int kappa = 0; kappa <= -2 + k; kappa++)
			{
				sumIkappa += -(vEG2G4G2c*(j - k + 2*n - 2*q - Max(0,-j + kappa))*(1 + j - k + Max(0,-j + kappa))*b[kappa])/2. - (vEG2G4G2cAWAH*(-1 - j + k - Max(0,-j + kappa))*(-j + k - 2*n + 2*q + Max(0,-j + kappa))*b[kappa]*psiexp(k - kappa,1))/2. - vEG2G4G2cAH*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-Max(0,-j + kappa),1)*((n - q)*psiexp(-j + k,1) + Max(0,-j + kappa)*(-1 + psiexp(1,1))*psiexp(-j + k,1) + (1 - n + q)*psiexp(1 - j + k,1) + (-1 - j + k - n + q)*psiexp(1 + Max(0,-j + kappa),1) + (j - k + n - q)*psiexp(2 + Max(0,-j + kappa),1)) + vEG2G4G2cAW*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(j - kappa,1)*((-j + k - n + q)*psiexp(-j + k,1) + (1 + j - k + n - q)*psiexp(1 - j + k,1) + (-1 + n - q)*psiexp(1 + Max(0,-j + kappa),1) + Max(0,-j + kappa)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,-j + kappa),1) + (-n + q)*psiexp(2 + Max(0,-j + kappa),1)) - vEG2G4G2cBH*b[kappa]*Power2(1/(-1 + psiexp(1,2)))*psiexp(-Max(0,-j + kappa),2)*((n - q)*psiexp(-j + k,2) + Max(0,-j + kappa)*(-1 + psiexp(1,2))*psiexp(-j + k,2) + (1 - n + q)*psiexp(1 - j + k,2) + (-1 - j + k - n + q)*psiexp(1 + Max(0,-j + kappa),2) + (j - k + n - q)*psiexp(2 + Max(0,-j + kappa),2)) - vEG2G4G2cAWBH*b[kappa]*Power2(1/(psiexp(1,1) - psiexp(1,2)))*psiexp(j - kappa,1)*psiexp(-Max(0,-j + kappa),2)*(-(Max(0,-j + kappa)*(psiexp(1,1) - psiexp(1,2))*psiexp(-j + k,2)*psiexp(1 + Max(0,-j + kappa),1)) + (1 - n + q)*psiexp(1 - j + k,2)*psiexp(1 + Max(0,-j + kappa),1) + (-1 - j + k - n + q)*psiexp(1 - j + k,1)*psiexp(1 + Max(0,-j + kappa),2) + (n - q)*psiexp(-j + k,2)*psiexp(2 + Max(0,-j + kappa),1) + (j - k + n - q)*psiexp(-j + k,1)*psiexp(2 + Max(0,-j + kappa),2)) - vEG2G4G2cAWCH*b[kappa]*Power2(1/(psiexp(1,1) - psiexp(1,3)))*psiexp(j - kappa,1)*psiexp(-Max(0,-j + kappa),3)*(-(Max(0,-j + kappa)*(psiexp(1,1) - psiexp(1,3))*psiexp(-j + k,3)*psiexp(1 + Max(0,-j + kappa),1)) + (1 - n + q)*psiexp(1 - j + k,3)*psiexp(1 + Max(0,-j + kappa),1) + (-1 - j + k - n + q)*psiexp(1 - j + k,1)*psiexp(1 + Max(0,-j + kappa),3) + (n - q)*psiexp(-j + k,3)*psiexp(2 + Max(0,-j + kappa),1) + (j - k + n - q)*psiexp(-j + k,1)*psiexp(2 + Max(0,-j + kappa),3));
			}
			MMjkP1 = sumIkappa;
		}
		MM(j, k) += (MMjkP1*b[j])/(n - q);
	}
}

for (int j = 1; j <= -2 + q; j++)
{
	for (int k = Max(1,2 + j); k <= q; k++)
	{
		MM(j, k) += -((1 + j - k)*(j - k + 2*n - 2*q)*vEG4G4c*b[j]*b[k])/(2.*(n - q)) + (vEG4G4cAD*b[j]*b[k]*Power2(1/(-1 + psiexp(1,1)))*(-n + q + (-1 + n - q)*psiexp(1,1) + (1 + j - k + n - q)*psiexp(1 + j - k,1) + (-j + k - n + q)*psiexp(2 + j - k,1))*psiexp(-j + k,1))/(n - q) + (vEG4G4cBD*b[j]*b[k]*Power2(1/(-1 + psiexp(1,2)))*(-n + q + (-1 + n - q)*psiexp(1,2) + (1 + j - k + n - q)*psiexp(1 + j - k,2) + (-j + k - n + q)*psiexp(2 + j - k,2))*psiexp(-j + k,2))/(n - q);
	}
}

for (int j = 1; j <= -2 + q; j++)
{
	for (int k = Max(1,2 + j); k <= q; k++)
	{
		scalar MMjkP1;
		{
			scalar sumIkappa = 0;
			for (int kappa = Max(0,1 + j); kappa <= -1 + k; kappa++)
			{
				sumIkappa += (j - kappa + n - q)*vEG6G2c*b[kappa] + (j - kappa + n - q)*vEG6G2cAH*b[kappa]*psiexp(k - kappa,1) + (j - kappa + n - q)*vEG6G2cBH*b[kappa]*psiexp(k - kappa,2) + (j - kappa + n - q)*vEG6G2cCH*b[kappa]*psiexp(k - kappa,3);
			}
			MMjkP1 = sumIkappa;
		}
		MM(j, k) += (MMjkP1*b[j])/(n - q);
	}
}

for (int j = 1; j <= -2 + q; j++)
{
	for (int k = 1; k <= -2 + q; k++)
	{
		scalar MMjkP1;
		{
			scalar sumIkappa = 0;
			for (int kappa = Max(2 + j,2 + k); kappa <= q; kappa++)
			{
				sumIkappa += -(vEG2G4G2c*(j - kappa + 2*n - 2*q - Max(0,-j + k))*(1 + j - kappa + Max(0,-j + k))*b[kappa])/2. - (vEG2G4G2cAWAH*(j - kappa + 2*n - 2*q - Max(0,-j + k))*(1 + j - kappa + Max(0,-j + k))*b[kappa]*psiexp(-k + kappa,1))/2. - vEG2G4G2cAH*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(-Max(0,-j + k),1)*((n - q)*psiexp(-j + kappa,1) + Max(0,-j + k)*(-1 + psiexp(1,1))*psiexp(-j + kappa,1) + (1 - n + q)*psiexp(1 - j + kappa,1) + (-1 - j + kappa - n + q)*psiexp(1 + Max(0,-j + k),1) + (j - kappa + n - q)*psiexp(2 + Max(0,-j + k),1)) + vEG2G4G2cAW*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(j - k,1)*((-j + kappa - n + q)*psiexp(-j + kappa,1) + (1 + j - kappa + n - q)*psiexp(1 - j + kappa,1) + (-1 + n - q)*psiexp(1 + Max(0,-j + k),1) + Max(0,-j + k)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,-j + k),1) + (-n + q)*psiexp(2 + Max(0,-j + k),1)) - vEG2G4G2cBH*b[kappa]*Power2(1/(-1 + psiexp(1,2)))*psiexp(-Max(0,-j + k),2)*((n - q)*psiexp(-j + kappa,2) + Max(0,-j + k)*(-1 + psiexp(1,2))*psiexp(-j + kappa,2) + (1 - n + q)*psiexp(1 - j + kappa,2) + (-1 - j + kappa - n + q)*psiexp(1 + Max(0,-j + k),2) + (j - kappa + n - q)*psiexp(2 + Max(0,-j + k),2)) - vEG2G4G2cAWBH*b[kappa]*Power2(1/(psiexp(1,1) - psiexp(1,2)))*psiexp(j - k,1)*psiexp(-Max(0,-j + k),2)*(-(Max(0,-j + k)*(psiexp(1,1) - psiexp(1,2))*psiexp(-j + kappa,2)*psiexp(1 + Max(0,-j + k),1)) + (1 - n + q)*psiexp(1 - j + kappa,2)*psiexp(1 + Max(0,-j + k),1) + (-1 - j + kappa - n + q)*psiexp(1 - j + kappa,1)*psiexp(1 + Max(0,-j + k),2) + (n - q)*psiexp(-j + kappa,2)*psiexp(2 + Max(0,-j + k),1) + (j - kappa + n - q)*psiexp(-j + kappa,1)*psiexp(2 + Max(0,-j + k),2)) - vEG2G4G2cAWCH*b[kappa]*Power2(1/(psiexp(1,1) - psiexp(1,3)))*psiexp(j - k,1)*psiexp(-Max(0,-j + k),3)*(-(Max(0,-j + k)*(psiexp(1,1) - psiexp(1,3))*psiexp(-j + kappa,3)*psiexp(1 + Max(0,-j + k),1)) + (1 - n + q)*psiexp(1 - j + kappa,3)*psiexp(1 + Max(0,-j + k),1) + (-1 - j + kappa - n + q)*psiexp(1 - j + kappa,1)*psiexp(1 + Max(0,-j + k),3) + (n - q)*psiexp(-j + kappa,3)*psiexp(2 + Max(0,-j + k),1) + (j - kappa + n - q)*psiexp(-j + kappa,1)*psiexp(2 + Max(0,-j + k),3));
			}
			MMjkP1 = sumIkappa;
		}
		MM(j, k) += (MMjkP1*b[j])/(n - q);
	}
}

for (int j = 1; j <= q; j++)
{
	for (int k = 1; k <= Min(-3 + j + n - q,-1 + q); k++)
	{
		scalar MMjkP1;
		{
			scalar sumIkappa = 0;
			for (int kappa = 1 + k; kappa <= Min(-2 + j + n - q,q); kappa++)
			{
				sumIkappa += (vEG2G2G4c*b[kappa]*(-n + q - 2*n*q + (1 - 2*n + 2*q)*Max(0,-j + kappa) + Power2(n) + Power2(q) + Power2(Max(0,-j + kappa))))/2. + (vEG2G2G4cAW*b[kappa]*(-n + q - 2*n*q + (1 - 2*n + 2*q)*Max(0,-j + kappa) + Power2(n) + Power2(q) + Power2(Max(0,-j + kappa)))*psiexp(-k + kappa,1))/2. + vEG2G2G4cAWAD*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(j - k,1)*(psiexp(1 + n - q,1) + (-1 + n - q)*psiexp(1 + Max(0,-j + kappa),1) + Max(0,-j + kappa)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,-j + kappa),1) + (-n + q)*psiexp(2 + Max(0,-j + kappa),1)) + vEG2G2G4cAD*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(j - kappa,1)*(psiexp(1 + n - q,1) + (-1 + n - q)*psiexp(1 + Max(0,-j + kappa),1) + Max(0,-j + kappa)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,-j + kappa),1) + (-n + q)*psiexp(2 + Max(0,-j + kappa),1)) + vEG2G2G4cAWBD*b[kappa]*Power2(1/(-1 + psiexp(1,2)))*psiexp(j - kappa,2)*psiexp(-k + kappa,1)*(psiexp(1 + n - q,2) + (-1 + n - q)*psiexp(1 + Max(0,-j + kappa),2) + Max(0,-j + kappa)*(-1 + psiexp(1,2))*psiexp(1 + Max(0,-j + kappa),2) + (-n + q)*psiexp(2 + Max(0,-j + kappa),2));
			}
			MMjkP1 = sumIkappa;
		}
		MM(j, k) += (MMjkP1*b[j])/(n - q);
	}
}

for (int j = 1; j <= q; j++)
{
	for (int k = 1; k <= Min(-2 + j + n - q,q); k++)
	{
		scalar MMjkP1;
		{
			scalar sumIkappa = 0;
			for (int kappa = 0; kappa <= -1 + k; kappa++)
			{
				sumIkappa += (vEG2G2G4c*b[kappa]*(-n + q - 2*n*q + (1 - 2*n + 2*q)*Max(0,-j + k) + Power2(n) + Power2(q) + Power2(Max(0,-j + k))))/2. + (vEG2G2G4cAW*b[kappa]*(-n + q - 2*n*q + (1 - 2*n + 2*q)*Max(0,-j + k) + Power2(n) + Power2(q) + Power2(Max(0,-j + k)))*psiexp(k - kappa,1))/2. + vEG2G2G4cAD*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(j - k,1)*(psiexp(1 + n - q,1) + (-1 + n - q)*psiexp(1 + Max(0,-j + k),1) + Max(0,-j + k)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,-j + k),1) + (-n + q)*psiexp(2 + Max(0,-j + k),1)) + vEG2G2G4cAWAD*b[kappa]*Power2(1/(-1 + psiexp(1,1)))*psiexp(j - kappa,1)*(psiexp(1 + n - q,1) + (-1 + n - q)*psiexp(1 + Max(0,-j + k),1) + Max(0,-j + k)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,-j + k),1) + (-n + q)*psiexp(2 + Max(0,-j + k),1)) + vEG2G2G4cAWBD*b[kappa]*Power2(1/(-1 + psiexp(1,2)))*psiexp(j - k,2)*psiexp(k - kappa,1)*(psiexp(1 + n - q,2) + (-1 + n - q)*psiexp(1 + Max(0,-j + k),2) + Max(0,-j + k)*(-1 + psiexp(1,2))*psiexp(1 + Max(0,-j + k),2) + (-n + q)*psiexp(2 + Max(0,-j + k),2));
			}
			MMjkP1 = sumIkappa;
		}
		MM(j, k) += (MMjkP1*b[j])/(n - q);
	}
}

for (int j = 1; j <= -1 + q; j++)
{
	for (int k = 1; k <= -1 + q; k++)
	{
		scalar MMjkP1;
		{
			scalar sumIkappa = 0;
			for (int kappa = Max(1 + j,1 + k); kappa <= q; kappa++)
			{
				sumIkappa += (j - kappa + n - q)*vEG2G6c*b[kappa] + (j - kappa + n - q)*vEG2G6cAW*b[kappa]*psiexp(-k + kappa,1);
			}
			MMjkP1 = sumIkappa;
		}
		MM(j, k) += (MMjkP1*b[j])/(n - q);
	}
}

for (int j = 1; j <= q; j++)
{
	for (int k = 1; k <= Min(-2 + j + n - q,q); k++)
	{
		MM(j, k) += (vEG4G4c*b[j]*b[k]*(-n + q - 2*n*q + (1 - 2*n + 2*q)*Max(0,-j + k) + Power2(n) + Power2(q) + Power2(Max(0,-j + k))))/(2.*(n - q)) + (vEG4G4cAD*b[j]*b[k]*Power2(1/(-1 + psiexp(1,1)))*psiexp(j - k,1)*(psiexp(1 + n - q,1) + (-1 + n - q)*psiexp(1 + Max(0,-j + k),1) + Max(0,-j + k)*(-1 + psiexp(1,1))*psiexp(1 + Max(0,-j + k),1) + (-n + q)*psiexp(2 + Max(0,-j + k),1)))/(n - q) + (vEG4G4cBD*b[j]*b[k]*Power2(1/(-1 + psiexp(1,2)))*psiexp(j - k,2)*(psiexp(1 + n - q,2) + (-1 + n - q)*psiexp(1 + Max(0,-j + k),2) + Max(0,-j + k)*(-1 + psiexp(1,2))*psiexp(1 + Max(0,-j + k),2) + (-n + q)*psiexp(2 + Max(0,-j + k),2)))/(n - q);
	}
}

for (int j = 1; j <= -2 + q; j++)
{
	for (int k = Max(1,1 + j); k <= -1 + q; k++)
	{
		scalar MMjkP1;
		{
			scalar sumIkappa = 0;
			for (int kappa = 1 + k; kappa <= q; kappa++)
			{
				sumIkappa += (j - k + n - q)*vEG6G2c*b[kappa] + (j - k + n - q)*vEG6G2cAH*b[kappa]*psiexp(-k + kappa,1) + (j - k + n - q)*vEG6G2cBH*b[kappa]*psiexp(-k + kappa,2) + (j - k + n - q)*vEG6G2cCH*b[kappa]*psiexp(-k + kappa,3);
			}
			MMjkP1 = sumIkappa;
		}
		MM(j, k) += (MMjkP1*b[j])/(n - q);
	}
}

for (int j = 1; j <= -1 + q; j++)
{
	for (int k = Max(1,1 + j); k <= q; k++)
	{
		scalar MMjkP1;
		{
			scalar sumIkappa = 0;
			for (int kappa = 0; kappa <= -1 + k; kappa++)
			{
				sumIkappa += (j - k + n - q)*vEG2G6c*b[kappa] + (j - k + n - q)*vEG2G6cAW*b[kappa]*psiexp(k - kappa,1);
			}
			MMjkP1 = sumIkappa;
		}
		MM(j, k) += (MMjkP1*b[j])/(n - q);
	}
}

for (int j = 1; j <= -1 + q; j++)
{
	for (int k = Max(1,1 + j); k <= q; k++)
	{
		MM(j, k) += ((j - k + n - q)*vEG8c*b[j]*b[k])/(n - q);
	}
}

free(psiexp_array);
#undef Psi
#undef psiexp
