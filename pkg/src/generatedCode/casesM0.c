#define psiexp(i, j) psiexp_array[3*(i - -2*q) + j - 1]
#define Psi(j) psiParam[j - 1]
const scalar h = *hParam;
scalar *psiexp_array = malloc(sizeof(scalar)*3*(1 + 4*q));

for (int i = -2*q; i <= 2*q; i++)
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
		{
			scalar sumInuIkappa = 0;
			for (int nu = 0; nu <= -1 + q; nu++)
			{
				for (int kappa = 1 + nu; kappa <= q; kappa++)
				{
					sumInuIkappa += vEG2G2c*b[kappa]*b[nu] + vEG2G2cAW*b[kappa]*b[nu]*psiexp(kappa - nu,1);
				}
			}
			MM(j, k) += sumInuIkappa;
		}
	}
}

for (int j = 0; j <= 0; j++)
{
	for (int k = 0; k <= 0; k++)
	{
		{
			scalar sumInuIkappa = 0;
			for (int nu = 1; nu <= q; nu++)
			{
				for (int kappa = 0; kappa <= -1 + nu; kappa++)
				{
					sumInuIkappa += vEG2G2c*b[kappa]*b[nu] + vEG2G2cAW*b[kappa]*b[nu]*psiexp(-kappa + nu,1);
				}
			}
			MM(j, k) += sumInuIkappa;
		}
	}
}

for (int j = 0; j <= 0; j++)
{
	for (int k = 0; k <= 0; k++)
	{
		{
			scalar sumInuIkappa = 0;
			for (int nu = 0; nu <= q; nu++)
			{
				for (int kappa = nu; kappa <= nu; kappa++)
				{
					sumInuIkappa += vEG4c*Power2(b[nu]);
				}
			}
			MM(j, k) += sumInuIkappa;
		}
	}
}

for (int j = 0; j <= 0; j++)
{
	for (int k = 0; k <= 0; k++)
	{
		{
			scalar sumIkappa = 0;
			for (int kappa = 0; kappa <= q; kappa++)
			{
				sumIkappa += 2*azero*vEG2c*b[kappa];
			}
			MM(j, k) += sumIkappa;
		}
	}
}

for (int j = 0; j <= 0; j++)
{
	for (int k = 0; k <= 0; k++)
	{
		MM(j, k) += Power2(azero);
	}
}

for (int j = 0; j <= 0; j++)
{
	for (int k = 1; k <= -2 + q; k++)
	{
		{
			scalar sumInuIkappa = 0;
			for (int nu = 1 + k; nu <= -1 + q; nu++)
			{
				for (int kappa = 1 + nu; kappa <= q; kappa++)
				{
					sumInuIkappa += vEG2G2G2c*b[kappa]*b[nu] + vEG2G2G2cAWAD*b[kappa]*b[nu]*psiexp(-k + kappa,1) + vEG2G2G2cAD*b[kappa]*b[nu]*psiexp(kappa - nu,1) + vEG2G2G2cAW*b[kappa]*b[nu]*psiexp(-k + nu,1) + vEG2G2G2cAWBD*b[kappa]*b[nu]*psiexp(kappa - nu,2)*psiexp(-k + nu,1);
				}
			}
			MM(j, k) += sumInuIkappa;
		}
	}
}

for (int j = 0; j <= 0; j++)
{
	for (int k = 1; k <= -2 + q; k++)
	{
		{
			scalar sumInuIkappa = 0;
			for (int nu = 2 + k; nu <= q; nu++)
			{
				for (int kappa = 1 + k; kappa <= -1 + nu; kappa++)
				{
					sumInuIkappa += vEG2G2G2c*b[kappa]*b[nu] + vEG2G2G2cAW*b[kappa]*b[nu]*psiexp(-k + kappa,1) + vEG2G2G2cAWAD*b[kappa]*b[nu]*psiexp(-k + nu,1) + vEG2G2G2cAD*b[kappa]*b[nu]*psiexp(-kappa + nu,1) + vEG2G2G2cAWBD*b[kappa]*b[nu]*psiexp(-k + kappa,1)*psiexp(-kappa + nu,2);
				}
			}
			MM(j, k) += sumInuIkappa;
		}
	}
}

for (int j = 0; j <= 0; j++)
{
	for (int k = 1; k <= -1 + q; k++)
	{
		{
			scalar sumInuIkappa = 0;
			for (int nu = 1 + k; nu <= q; nu++)
			{
				for (int kappa = 0; kappa <= -1 + k; kappa++)
				{
					sumInuIkappa += vEG2G2G2c*b[kappa]*b[nu] + vEG2G2G2cAW*b[kappa]*b[nu]*psiexp(k - kappa,1) + vEG2G2G2cAD*b[kappa]*b[nu]*psiexp(-k + nu,1) + vEG2G2G2cAWBD*b[kappa]*b[nu]*psiexp(k - kappa,1)*psiexp(-k + nu,2) + vEG2G2G2cAWAD*b[kappa]*b[nu]*psiexp(-kappa + nu,1);
				}
			}
			MM(j, k) += sumInuIkappa;
		}
	}
}

for (int j = 0; j <= 0; j++)
{
	for (int k = 1; k <= -1 + q; k++)
	{
		{
			scalar sumInuIkappa = 0;
			for (int nu = 1 + k; nu <= q; nu++)
			{
				for (int kappa = nu; kappa <= nu; kappa++)
				{
					sumInuIkappa += vEG2G4c*Power2(b[nu]) + vEG2G4cAW*Power2(b[nu])*psiexp(-k + nu,1);
				}
			}
			MM(j, k) += sumInuIkappa;
		}
	}
}

for (int j = 0; j <= 0; j++)
{
	for (int k = 1; k <= -1 + q; k++)
	{
		{
			scalar sumInuIkappa = 0;
			for (int nu = 1 + k; nu <= q; nu++)
			{
				for (int kappa = k; kappa <= k; kappa++)
				{
					sumInuIkappa += vEG4G2c*b[kappa]*b[nu] + vEG4G2cAD*b[kappa]*b[nu]*psiexp(-kappa + nu,1) + vEG4G2cBD*b[kappa]*b[nu]*psiexp(-kappa + nu,2);
				}
			}
			MM(j, k) += sumInuIkappa;
		}
	}
}

for (int j = 0; j <= 0; j++)
{
	for (int k = 1; k <= -1 + q; k++)
	{
		{
			scalar sumInuIkappa = 0;
			for (int nu = 0; nu <= -1 + k; nu++)
			{
				for (int kappa = 1 + k; kappa <= q; kappa++)
				{
					sumInuIkappa += vEG2G2G2c*b[kappa]*b[nu] + vEG2G2G2cAD*b[kappa]*b[nu]*psiexp(-k + kappa,1) + vEG2G2G2cAW*b[kappa]*b[nu]*psiexp(k - nu,1) + vEG2G2G2cAWBD*b[kappa]*b[nu]*psiexp(-k + kappa,2)*psiexp(k - nu,1) + vEG2G2G2cAWAD*b[kappa]*b[nu]*psiexp(kappa - nu,1);
				}
			}
			MM(j, k) += sumInuIkappa;
		}
	}
}

for (int j = 0; j <= 0; j++)
{
	for (int k = 2; k <= q; k++)
	{
		{
			scalar sumInuIkappa = 0;
			for (int nu = 0; nu <= -2 + k; nu++)
			{
				for (int kappa = 1 + nu; kappa <= -1 + k; kappa++)
				{
					sumInuIkappa += vEG2G2G2c*b[kappa]*b[nu] + vEG2G2G2cAD*b[kappa]*b[nu]*psiexp(k - kappa,1) + vEG2G2G2cAWAD*b[kappa]*b[nu]*psiexp(k - nu,1) + vEG2G2G2cAW*b[kappa]*b[nu]*psiexp(kappa - nu,1) + vEG2G2G2cAWBD*b[kappa]*b[nu]*psiexp(k - kappa,2)*psiexp(kappa - nu,1);
				}
			}
			MM(j, k) += sumInuIkappa;
		}
	}
}

for (int j = 0; j <= 0; j++)
{
	for (int k = 2; k <= q; k++)
	{
		{
			scalar sumInuIkappa = 0;
			for (int nu = 1; nu <= -1 + k; nu++)
			{
				for (int kappa = 0; kappa <= -1 + nu; kappa++)
				{
					sumInuIkappa += vEG2G2G2c*b[kappa]*b[nu] + vEG2G2G2cAWAD*b[kappa]*b[nu]*psiexp(k - kappa,1) + vEG2G2G2cAD*b[kappa]*b[nu]*psiexp(k - nu,1) + vEG2G2G2cAW*b[kappa]*b[nu]*psiexp(-kappa + nu,1) + vEG2G2G2cAWBD*b[kappa]*b[nu]*psiexp(k - nu,2)*psiexp(-kappa + nu,1);
				}
			}
			MM(j, k) += sumInuIkappa;
		}
	}
}

for (int j = 0; j <= 0; j++)
{
	for (int k = 1; k <= q; k++)
	{
		{
			scalar sumInuIkappa = 0;
			for (int nu = 0; nu <= -1 + k; nu++)
			{
				for (int kappa = k; kappa <= k; kappa++)
				{
					sumInuIkappa += vEG2G4c*b[kappa]*b[nu] + vEG2G4cAW*b[kappa]*b[nu]*psiexp(kappa - nu,1);
				}
			}
			MM(j, k) += sumInuIkappa;
		}
	}
}

for (int j = 0; j <= 0; j++)
{
	for (int k = 1; k <= q; k++)
	{
		{
			scalar sumInuIkappa = 0;
			for (int nu = 0; nu <= -1 + k; nu++)
			{
				for (int kappa = nu; kappa <= nu; kappa++)
				{
					sumInuIkappa += vEG4G2c*Power2(b[nu]) + vEG4G2cAD*Power2(b[nu])*psiexp(k - nu,1) + vEG4G2cBD*Power2(b[nu])*psiexp(k - nu,2);
				}
			}
			MM(j, k) += sumInuIkappa;
		}
	}
}

for (int j = 0; j <= 0; j++)
{
	for (int k = 1; k <= -1 + q; k++)
	{
		{
			scalar sumInuIkappa = 0;
			for (int nu = k; nu <= k; nu++)
			{
				for (int kappa = 1 + k; kappa <= q; kappa++)
				{
					sumInuIkappa += vEG4G2c*b[kappa]*b[nu] + vEG4G2cAD*b[kappa]*b[nu]*psiexp(kappa - nu,1) + vEG4G2cBD*b[kappa]*b[nu]*psiexp(kappa - nu,2);
				}
			}
			MM(j, k) += sumInuIkappa;
		}
	}
}

for (int j = 0; j <= 0; j++)
{
	for (int k = 1; k <= q; k++)
	{
		{
			scalar sumInuIkappa = 0;
			for (int nu = k; nu <= k; nu++)
			{
				for (int kappa = 0; kappa <= -1 + k; kappa++)
				{
					sumInuIkappa += vEG2G4c*b[kappa]*b[nu] + vEG2G4cAW*b[kappa]*b[nu]*psiexp(-kappa + nu,1);
				}
			}
			MM(j, k) += sumInuIkappa;
		}
	}
}

for (int j = 0; j <= 0; j++)
{
	for (int k = 1; k <= q; k++)
	{
		{
			scalar sumInuIkappa = 0;
			for (int nu = k; nu <= k; nu++)
			{
				for (int kappa = k; kappa <= k; kappa++)
				{
					sumInuIkappa += vEG6c*Power2(b[nu]);
				}
			}
			MM(j, k) += sumInuIkappa;
		}
	}
}

for (int j = 0; j <= 0; j++)
{
	for (int k = 1; k <= -1 + q; k++)
	{
		{
			scalar sumInu = 0;
			for (int nu = 1 + k; nu <= q; nu++)
			{
				sumInu += azero*vEG2G2c*b[nu] + azero*vEG2G2cAW*b[nu]*psiexp(-k + nu,1);
			}
			MM(j, k) += sumInu;
		}
	}
}

for (int j = 0; j <= 0; j++)
{
	for (int k = 1; k <= q; k++)
	{
		{
			scalar sumInu = 0;
			for (int nu = 0; nu <= -1 + k; nu++)
			{
				sumInu += azero*vEG2G2c*b[nu] + azero*vEG2G2cAW*b[nu]*psiexp(k - nu,1);
			}
			MM(j, k) += sumInu;
		}
	}
}

for (int j = 0; j <= 0; j++)
{
	for (int k = 1; k <= q; k++)
	{
		{
			scalar sumInu = 0;
			for (int nu = k; nu <= k; nu++)
			{
				sumInu += azero*vEG4c*b[nu];
			}
			MM(j, k) += sumInu;
		}
	}
}

for (int j = 0; j <= 0; j++)
{
	for (int k = 1; k <= -1 + q; k++)
	{
		{
			scalar sumInu = 0;
			for (int nu = 1 + k; nu <= q; nu++)
			{
				sumInu += azero*vEG2G2c*b[nu] + azero*vEG2G2cAW*b[nu]*psiexp(-k + nu,1);
			}
			MM(j, k) += sumInu;
		}
	}
}

for (int j = 0; j <= 0; j++)
{
	for (int k = 1; k <= q; k++)
	{
		{
			scalar sumInu = 0;
			for (int nu = 0; nu <= -1 + k; nu++)
			{
				sumInu += azero*vEG2G2c*b[nu] + azero*vEG2G2cAW*b[nu]*psiexp(k - nu,1);
			}
			MM(j, k) += sumInu;
		}
	}
}

for (int j = 0; j <= 0; j++)
{
	for (int k = 1; k <= q; k++)
	{
		{
			scalar sumInu = 0;
			for (int nu = k; nu <= k; nu++)
			{
				sumInu += azero*vEG4c*b[nu];
			}
			MM(j, k) += sumInu;
		}
	}
}

for (int j = 0; j <= 0; j++)
{
	for (int k = 1; k <= q; k++)
	{
		MM(j, k) += vEG2c*Power2(azero);
	}
}

for (int j = 1; j <= -2 + q; j++)
{
	for (int k = 0; k <= 0; k++)
	{
		{
			scalar sumInuIkappa = 0;
			for (int nu = 1 + j; nu <= -1 + q; nu++)
			{
				for (int kappa = 1 + nu; kappa <= q; kappa++)
				{
					sumInuIkappa += vEG2G2G2c*b[kappa]*b[nu] + vEG2G2G2cAWAD*b[kappa]*b[nu]*psiexp(-j + kappa,1) + vEG2G2G2cAD*b[kappa]*b[nu]*psiexp(kappa - nu,1) + vEG2G2G2cAW*b[kappa]*b[nu]*psiexp(-j + nu,1) + vEG2G2G2cAWBD*b[kappa]*b[nu]*psiexp(kappa - nu,2)*psiexp(-j + nu,1);
				}
			}
			MM(j, k) += sumInuIkappa;
		}
	}
}

for (int j = 1; j <= -2 + q; j++)
{
	for (int k = 0; k <= 0; k++)
	{
		{
			scalar sumInuIkappa = 0;
			for (int nu = 2 + j; nu <= q; nu++)
			{
				for (int kappa = 1 + j; kappa <= -1 + nu; kappa++)
				{
					sumInuIkappa += vEG2G2G2c*b[kappa]*b[nu] + vEG2G2G2cAW*b[kappa]*b[nu]*psiexp(-j + kappa,1) + vEG2G2G2cAWAD*b[kappa]*b[nu]*psiexp(-j + nu,1) + vEG2G2G2cAD*b[kappa]*b[nu]*psiexp(-kappa + nu,1) + vEG2G2G2cAWBD*b[kappa]*b[nu]*psiexp(-j + kappa,1)*psiexp(-kappa + nu,2);
				}
			}
			MM(j, k) += sumInuIkappa;
		}
	}
}

for (int j = 1; j <= -1 + q; j++)
{
	for (int k = 0; k <= 0; k++)
	{
		{
			scalar sumInuIkappa = 0;
			for (int nu = 1 + j; nu <= q; nu++)
			{
				for (int kappa = 0; kappa <= -1 + j; kappa++)
				{
					sumInuIkappa += vEG2G2G2c*b[kappa]*b[nu] + vEG2G2G2cAW*b[kappa]*b[nu]*psiexp(j - kappa,1) + vEG2G2G2cAD*b[kappa]*b[nu]*psiexp(-j + nu,1) + vEG2G2G2cAWBD*b[kappa]*b[nu]*psiexp(j - kappa,1)*psiexp(-j + nu,2) + vEG2G2G2cAWAD*b[kappa]*b[nu]*psiexp(-kappa + nu,1);
				}
			}
			MM(j, k) += sumInuIkappa;
		}
	}
}

for (int j = 1; j <= -1 + q; j++)
{
	for (int k = 0; k <= 0; k++)
	{
		{
			scalar sumInuIkappa = 0;
			for (int nu = 1 + j; nu <= q; nu++)
			{
				for (int kappa = nu; kappa <= nu; kappa++)
				{
					sumInuIkappa += vEG2G4c*Power2(b[nu]) + vEG2G4cAW*Power2(b[nu])*psiexp(-j + nu,1);
				}
			}
			MM(j, k) += sumInuIkappa;
		}
	}
}

for (int j = 1; j <= -1 + q; j++)
{
	for (int k = 0; k <= 0; k++)
	{
		{
			scalar sumInuIkappa = 0;
			for (int nu = 1 + j; nu <= q; nu++)
			{
				for (int kappa = j; kappa <= j; kappa++)
				{
					sumInuIkappa += vEG4G2c*b[kappa]*b[nu] + vEG4G2cAD*b[kappa]*b[nu]*psiexp(-kappa + nu,1) + vEG4G2cBD*b[kappa]*b[nu]*psiexp(-kappa + nu,2);
				}
			}
			MM(j, k) += sumInuIkappa;
		}
	}
}

for (int j = 1; j <= -1 + q; j++)
{
	for (int k = 0; k <= 0; k++)
	{
		{
			scalar sumInuIkappa = 0;
			for (int nu = 0; nu <= -1 + j; nu++)
			{
				for (int kappa = 1 + j; kappa <= q; kappa++)
				{
					sumInuIkappa += vEG2G2G2c*b[kappa]*b[nu] + vEG2G2G2cAD*b[kappa]*b[nu]*psiexp(-j + kappa,1) + vEG2G2G2cAW*b[kappa]*b[nu]*psiexp(j - nu,1) + vEG2G2G2cAWBD*b[kappa]*b[nu]*psiexp(-j + kappa,2)*psiexp(j - nu,1) + vEG2G2G2cAWAD*b[kappa]*b[nu]*psiexp(kappa - nu,1);
				}
			}
			MM(j, k) += sumInuIkappa;
		}
	}
}

for (int j = 2; j <= q; j++)
{
	for (int k = 0; k <= 0; k++)
	{
		{
			scalar sumInuIkappa = 0;
			for (int nu = 0; nu <= -2 + j; nu++)
			{
				for (int kappa = 1 + nu; kappa <= -1 + j; kappa++)
				{
					sumInuIkappa += vEG2G2G2c*b[kappa]*b[nu] + vEG2G2G2cAD*b[kappa]*b[nu]*psiexp(j - kappa,1) + vEG2G2G2cAWAD*b[kappa]*b[nu]*psiexp(j - nu,1) + vEG2G2G2cAW*b[kappa]*b[nu]*psiexp(kappa - nu,1) + vEG2G2G2cAWBD*b[kappa]*b[nu]*psiexp(j - kappa,2)*psiexp(kappa - nu,1);
				}
			}
			MM(j, k) += sumInuIkappa;
		}
	}
}

for (int j = 2; j <= q; j++)
{
	for (int k = 0; k <= 0; k++)
	{
		{
			scalar sumInuIkappa = 0;
			for (int nu = 1; nu <= -1 + j; nu++)
			{
				for (int kappa = 0; kappa <= -1 + nu; kappa++)
				{
					sumInuIkappa += vEG2G2G2c*b[kappa]*b[nu] + vEG2G2G2cAWAD*b[kappa]*b[nu]*psiexp(j - kappa,1) + vEG2G2G2cAD*b[kappa]*b[nu]*psiexp(j - nu,1) + vEG2G2G2cAW*b[kappa]*b[nu]*psiexp(-kappa + nu,1) + vEG2G2G2cAWBD*b[kappa]*b[nu]*psiexp(j - nu,2)*psiexp(-kappa + nu,1);
				}
			}
			MM(j, k) += sumInuIkappa;
		}
	}
}

for (int j = 1; j <= q; j++)
{
	for (int k = 0; k <= 0; k++)
	{
		{
			scalar sumInuIkappa = 0;
			for (int nu = 0; nu <= -1 + j; nu++)
			{
				for (int kappa = j; kappa <= j; kappa++)
				{
					sumInuIkappa += vEG2G4c*b[kappa]*b[nu] + vEG2G4cAW*b[kappa]*b[nu]*psiexp(kappa - nu,1);
				}
			}
			MM(j, k) += sumInuIkappa;
		}
	}
}

for (int j = 1; j <= q; j++)
{
	for (int k = 0; k <= 0; k++)
	{
		{
			scalar sumInuIkappa = 0;
			for (int nu = 0; nu <= -1 + j; nu++)
			{
				for (int kappa = nu; kappa <= nu; kappa++)
				{
					sumInuIkappa += vEG4G2c*Power2(b[nu]) + vEG4G2cAD*Power2(b[nu])*psiexp(j - nu,1) + vEG4G2cBD*Power2(b[nu])*psiexp(j - nu,2);
				}
			}
			MM(j, k) += sumInuIkappa;
		}
	}
}

for (int j = 1; j <= -1 + q; j++)
{
	for (int k = 0; k <= 0; k++)
	{
		{
			scalar sumInuIkappa = 0;
			for (int nu = j; nu <= j; nu++)
			{
				for (int kappa = 1 + j; kappa <= q; kappa++)
				{
					sumInuIkappa += vEG4G2c*b[kappa]*b[nu] + vEG4G2cAD*b[kappa]*b[nu]*psiexp(kappa - nu,1) + vEG4G2cBD*b[kappa]*b[nu]*psiexp(kappa - nu,2);
				}
			}
			MM(j, k) += sumInuIkappa;
		}
	}
}

for (int j = 1; j <= q; j++)
{
	for (int k = 0; k <= 0; k++)
	{
		{
			scalar sumInuIkappa = 0;
			for (int nu = j; nu <= j; nu++)
			{
				for (int kappa = 0; kappa <= -1 + j; kappa++)
				{
					sumInuIkappa += vEG2G4c*b[kappa]*b[nu] + vEG2G4cAW*b[kappa]*b[nu]*psiexp(-kappa + nu,1);
				}
			}
			MM(j, k) += sumInuIkappa;
		}
	}
}

for (int j = 1; j <= q; j++)
{
	for (int k = 0; k <= 0; k++)
	{
		{
			scalar sumInuIkappa = 0;
			for (int nu = j; nu <= j; nu++)
			{
				for (int kappa = j; kappa <= j; kappa++)
				{
					sumInuIkappa += vEG6c*Power2(b[nu]);
				}
			}
			MM(j, k) += sumInuIkappa;
		}
	}
}

for (int j = 1; j <= -1 + q; j++)
{
	for (int k = 0; k <= 0; k++)
	{
		{
			scalar sumIkappa = 0;
			for (int kappa = 1 + j; kappa <= q; kappa++)
			{
				sumIkappa += azero*vEG2G2c*b[kappa] + azero*vEG2G2cAW*b[kappa]*psiexp(-j + kappa,1);
			}
			MM(j, k) += sumIkappa;
		}
	}
}

for (int j = 1; j <= q; j++)
{
	for (int k = 0; k <= 0; k++)
	{
		{
			scalar sumIkappa = 0;
			for (int kappa = 0; kappa <= -1 + j; kappa++)
			{
				sumIkappa += azero*vEG2G2c*b[kappa] + azero*vEG2G2cAW*b[kappa]*psiexp(j - kappa,1);
			}
			MM(j, k) += sumIkappa;
		}
	}
}

for (int j = 1; j <= q; j++)
{
	for (int k = 0; k <= 0; k++)
	{
		{
			scalar sumIkappa = 0;
			for (int kappa = j; kappa <= j; kappa++)
			{
				sumIkappa += azero*vEG4c*b[kappa];
			}
			MM(j, k) += sumIkappa;
		}
	}
}

for (int j = 1; j <= -1 + q; j++)
{
	for (int k = 0; k <= 0; k++)
	{
		{
			scalar sumIkappa = 0;
			for (int kappa = 1 + j; kappa <= q; kappa++)
			{
				sumIkappa += azero*vEG2G2c*b[kappa] + azero*vEG2G2cAW*b[kappa]*psiexp(-j + kappa,1);
			}
			MM(j, k) += sumIkappa;
		}
	}
}

for (int j = 1; j <= q; j++)
{
	for (int k = 0; k <= 0; k++)
	{
		{
			scalar sumIkappa = 0;
			for (int kappa = 0; kappa <= -1 + j; kappa++)
			{
				sumIkappa += azero*vEG2G2c*b[kappa] + azero*vEG2G2cAW*b[kappa]*psiexp(j - kappa,1);
			}
			MM(j, k) += sumIkappa;
		}
	}
}

for (int j = 1; j <= q; j++)
{
	for (int k = 0; k <= 0; k++)
	{
		{
			scalar sumIkappa = 0;
			for (int kappa = j; kappa <= j; kappa++)
			{
				sumIkappa += azero*vEG4c*b[kappa];
			}
			MM(j, k) += sumIkappa;
		}
	}
}

for (int j = 1; j <= q; j++)
{
	for (int k = 0; k <= 0; k++)
	{
		MM(j, k) += vEG2c*Power2(azero);
	}
}

for (int j = 1; j <= -1 + q; j++)
{
	for (int k = 1 + j; k <= q; k++)
	{
		MM(j, k) += vEG2G2c*Power2(azero) + vEG2G2cAW*Power2(azero)*psiexp(-j + k,1);
	}
}

for (int j = 2; j <= q; j++)
{
	for (int k = 1; k <= -1 + j; k++)
	{
		MM(j, k) += vEG2G2c*Power2(azero) + vEG2G2cAW*Power2(azero)*psiexp(j - k,1);
	}
}

for (int j = 1; j <= q; j++)
{
	for (int k = j; k <= j; k++)
	{
		MM(j, k) += vEG4c*Power2(azero);
	}
}

for (int j = 1; j <= -2 + q; j++)
{
	for (int k = 1 + j; k <= -1 + q; k++)
	{
		{
			scalar sumInu = 0;
			for (int nu = 1 + k; nu <= q; nu++)
			{
				sumInu += azero*vEG2G2G2c*b[nu] + azero*vEG2G2G2cAW*b[nu]*psiexp(-j + k,1) + azero*vEG2G2G2cAWAD*b[nu]*psiexp(-j + nu,1) + azero*vEG2G2G2cAD*b[nu]*psiexp(-k + nu,1) + azero*vEG2G2G2cAWBD*b[nu]*psiexp(-j + k,1)*psiexp(-k + nu,2);
			}
			MM(j, k) += sumInu;
		}
	}
}

for (int j = 1; j <= -2 + q; j++)
{
	for (int k = 2 + j; k <= q; k++)
	{
		{
			scalar sumInu = 0;
			for (int nu = 1 + j; nu <= -1 + k; nu++)
			{
				sumInu += azero*vEG2G2G2c*b[nu] + azero*vEG2G2G2cAWAD*b[nu]*psiexp(-j + k,1) + azero*vEG2G2G2cAD*b[nu]*psiexp(k - nu,1) + azero*vEG2G2G2cAW*b[nu]*psiexp(-j + nu,1) + azero*vEG2G2G2cAWBD*b[nu]*psiexp(k - nu,2)*psiexp(-j + nu,1);
			}
			MM(j, k) += sumInu;
		}
	}
}

for (int j = 1; j <= -1 + q; j++)
{
	for (int k = 1 + j; k <= q; k++)
	{
		{
			scalar sumInu = 0;
			for (int nu = 0; nu <= -1 + j; nu++)
			{
				sumInu += azero*vEG2G2G2c*b[nu] + azero*vEG2G2G2cAD*b[nu]*psiexp(-j + k,1) + azero*vEG2G2G2cAW*b[nu]*psiexp(j - nu,1) + azero*vEG2G2G2cAWBD*b[nu]*psiexp(-j + k,2)*psiexp(j - nu,1) + azero*vEG2G2G2cAWAD*b[nu]*psiexp(k - nu,1);
			}
			MM(j, k) += sumInu;
		}
	}
}

for (int j = 1; j <= -1 + q; j++)
{
	for (int k = 1 + j; k <= q; k++)
	{
		{
			scalar sumInu = 0;
			for (int nu = k; nu <= k; nu++)
			{
				sumInu += azero*vEG2G4c*b[nu] + azero*vEG2G4cAW*b[nu]*psiexp(-j + nu,1);
			}
			MM(j, k) += sumInu;
		}
	}
}

for (int j = 1; j <= -1 + q; j++)
{
	for (int k = 1 + j; k <= q; k++)
	{
		{
			scalar sumInu = 0;
			for (int nu = j; nu <= j; nu++)
			{
				sumInu += azero*vEG4G2c*b[nu] + azero*vEG4G2cAD*b[nu]*psiexp(k - nu,1) + azero*vEG4G2cBD*b[nu]*psiexp(k - nu,2);
			}
			MM(j, k) += sumInu;
		}
	}
}

for (int j = 2; j <= -1 + q; j++)
{
	for (int k = 1; k <= -1 + j; k++)
	{
		{
			scalar sumInu = 0;
			for (int nu = 1 + j; nu <= q; nu++)
			{
				sumInu += azero*vEG2G2G2c*b[nu] + azero*vEG2G2G2cAW*b[nu]*psiexp(j - k,1) + azero*vEG2G2G2cAD*b[nu]*psiexp(-j + nu,1) + azero*vEG2G2G2cAWBD*b[nu]*psiexp(j - k,1)*psiexp(-j + nu,2) + azero*vEG2G2G2cAWAD*b[nu]*psiexp(-k + nu,1);
			}
			MM(j, k) += sumInu;
		}
	}
}

for (int j = 3; j <= q; j++)
{
	for (int k = 1; k <= -2 + j; k++)
	{
		{
			scalar sumInu = 0;
			for (int nu = 1 + k; nu <= -1 + j; nu++)
			{
				sumInu += azero*vEG2G2G2c*b[nu] + azero*vEG2G2G2cAWAD*b[nu]*psiexp(j - k,1) + azero*vEG2G2G2cAD*b[nu]*psiexp(j - nu,1) + azero*vEG2G2G2cAW*b[nu]*psiexp(-k + nu,1) + azero*vEG2G2G2cAWBD*b[nu]*psiexp(j - nu,2)*psiexp(-k + nu,1);
			}
			MM(j, k) += sumInu;
		}
	}
}

for (int j = 2; j <= q; j++)
{
	for (int k = 1; k <= -1 + j; k++)
	{
		{
			scalar sumInu = 0;
			for (int nu = 0; nu <= -1 + k; nu++)
			{
				sumInu += azero*vEG2G2G2c*b[nu] + azero*vEG2G2G2cAD*b[nu]*psiexp(j - k,1) + azero*vEG2G2G2cAWAD*b[nu]*psiexp(j - nu,1) + azero*vEG2G2G2cAW*b[nu]*psiexp(k - nu,1) + azero*vEG2G2G2cAWBD*b[nu]*psiexp(j - k,2)*psiexp(k - nu,1);
			}
			MM(j, k) += sumInu;
		}
	}
}

for (int j = 2; j <= q; j++)
{
	for (int k = 1; k <= -1 + j; k++)
	{
		{
			scalar sumInu = 0;
			for (int nu = j; nu <= j; nu++)
			{
				sumInu += azero*vEG2G4c*b[nu] + azero*vEG2G4cAW*b[nu]*psiexp(-k + nu,1);
			}
			MM(j, k) += sumInu;
		}
	}
}

for (int j = 2; j <= q; j++)
{
	for (int k = 1; k <= -1 + j; k++)
	{
		{
			scalar sumInu = 0;
			for (int nu = k; nu <= k; nu++)
			{
				sumInu += azero*vEG4G2c*b[nu] + azero*vEG4G2cAD*b[nu]*psiexp(j - nu,1) + azero*vEG4G2cBD*b[nu]*psiexp(j - nu,2);
			}
			MM(j, k) += sumInu;
		}
	}
}

for (int j = 1; j <= -1 + q; j++)
{
	for (int k = j; k <= j; k++)
	{
		{
			scalar sumInu = 0;
			for (int nu = 1 + j; nu <= q; nu++)
			{
				sumInu += azero*vEG4G2c*b[nu] + azero*vEG4G2cAD*b[nu]*psiexp(-k + nu,1) + azero*vEG4G2cBD*b[nu]*psiexp(-k + nu,2);
			}
			MM(j, k) += sumInu;
		}
	}
}

for (int j = 1; j <= q; j++)
{
	for (int k = j; k <= j; k++)
	{
		{
			scalar sumInu = 0;
			for (int nu = 0; nu <= -1 + j; nu++)
			{
				sumInu += azero*vEG2G4c*b[nu] + azero*vEG2G4cAW*b[nu]*psiexp(k - nu,1);
			}
			MM(j, k) += sumInu;
		}
	}
}

for (int j = 1; j <= q; j++)
{
	for (int k = j; k <= j; k++)
	{
		{
			scalar sumInu = 0;
			for (int nu = j; nu <= j; nu++)
			{
				sumInu += azero*vEG6c*b[nu];
			}
			MM(j, k) += sumInu;
		}
	}
}

for (int j = 1; j <= -2 + q; j++)
{
	for (int k = 1 + j; k <= -1 + q; k++)
	{
		{
			scalar sumInu = 0;
			for (int nu = 1 + k; nu <= q; nu++)
			{
				sumInu += azero*vEG2G2G2c*b[nu] + azero*vEG2G2G2cAW*b[nu]*psiexp(-j + k,1) + azero*vEG2G2G2cAWAD*b[nu]*psiexp(-j + nu,1) + azero*vEG2G2G2cAD*b[nu]*psiexp(-k + nu,1) + azero*vEG2G2G2cAWBD*b[nu]*psiexp(-j + k,1)*psiexp(-k + nu,2);
			}
			MM(j, k) += sumInu;
		}
	}
}

for (int j = 1; j <= -2 + q; j++)
{
	for (int k = 2 + j; k <= q; k++)
	{
		{
			scalar sumInu = 0;
			for (int nu = 1 + j; nu <= -1 + k; nu++)
			{
				sumInu += azero*vEG2G2G2c*b[nu] + azero*vEG2G2G2cAWAD*b[nu]*psiexp(-j + k,1) + azero*vEG2G2G2cAD*b[nu]*psiexp(k - nu,1) + azero*vEG2G2G2cAW*b[nu]*psiexp(-j + nu,1) + azero*vEG2G2G2cAWBD*b[nu]*psiexp(k - nu,2)*psiexp(-j + nu,1);
			}
			MM(j, k) += sumInu;
		}
	}
}

for (int j = 1; j <= -1 + q; j++)
{
	for (int k = 1 + j; k <= q; k++)
	{
		{
			scalar sumInu = 0;
			for (int nu = 0; nu <= -1 + j; nu++)
			{
				sumInu += azero*vEG2G2G2c*b[nu] + azero*vEG2G2G2cAD*b[nu]*psiexp(-j + k,1) + azero*vEG2G2G2cAW*b[nu]*psiexp(j - nu,1) + azero*vEG2G2G2cAWBD*b[nu]*psiexp(-j + k,2)*psiexp(j - nu,1) + azero*vEG2G2G2cAWAD*b[nu]*psiexp(k - nu,1);
			}
			MM(j, k) += sumInu;
		}
	}
}

for (int j = 1; j <= -1 + q; j++)
{
	for (int k = 1 + j; k <= q; k++)
	{
		{
			scalar sumInu = 0;
			for (int nu = k; nu <= k; nu++)
			{
				sumInu += azero*vEG2G4c*b[nu] + azero*vEG2G4cAW*b[nu]*psiexp(-j + nu,1);
			}
			MM(j, k) += sumInu;
		}
	}
}

for (int j = 1; j <= -1 + q; j++)
{
	for (int k = 1 + j; k <= q; k++)
	{
		{
			scalar sumInu = 0;
			for (int nu = j; nu <= j; nu++)
			{
				sumInu += azero*vEG4G2c*b[nu] + azero*vEG4G2cAD*b[nu]*psiexp(k - nu,1) + azero*vEG4G2cBD*b[nu]*psiexp(k - nu,2);
			}
			MM(j, k) += sumInu;
		}
	}
}

for (int j = 2; j <= -1 + q; j++)
{
	for (int k = 1; k <= -1 + j; k++)
	{
		{
			scalar sumInu = 0;
			for (int nu = 1 + j; nu <= q; nu++)
			{
				sumInu += azero*vEG2G2G2c*b[nu] + azero*vEG2G2G2cAW*b[nu]*psiexp(j - k,1) + azero*vEG2G2G2cAD*b[nu]*psiexp(-j + nu,1) + azero*vEG2G2G2cAWBD*b[nu]*psiexp(j - k,1)*psiexp(-j + nu,2) + azero*vEG2G2G2cAWAD*b[nu]*psiexp(-k + nu,1);
			}
			MM(j, k) += sumInu;
		}
	}
}

for (int j = 3; j <= q; j++)
{
	for (int k = 1; k <= -2 + j; k++)
	{
		{
			scalar sumInu = 0;
			for (int nu = 1 + k; nu <= -1 + j; nu++)
			{
				sumInu += azero*vEG2G2G2c*b[nu] + azero*vEG2G2G2cAWAD*b[nu]*psiexp(j - k,1) + azero*vEG2G2G2cAD*b[nu]*psiexp(j - nu,1) + azero*vEG2G2G2cAW*b[nu]*psiexp(-k + nu,1) + azero*vEG2G2G2cAWBD*b[nu]*psiexp(j - nu,2)*psiexp(-k + nu,1);
			}
			MM(j, k) += sumInu;
		}
	}
}

for (int j = 2; j <= q; j++)
{
	for (int k = 1; k <= -1 + j; k++)
	{
		{
			scalar sumInu = 0;
			for (int nu = 0; nu <= -1 + k; nu++)
			{
				sumInu += azero*vEG2G2G2c*b[nu] + azero*vEG2G2G2cAD*b[nu]*psiexp(j - k,1) + azero*vEG2G2G2cAWAD*b[nu]*psiexp(j - nu,1) + azero*vEG2G2G2cAW*b[nu]*psiexp(k - nu,1) + azero*vEG2G2G2cAWBD*b[nu]*psiexp(j - k,2)*psiexp(k - nu,1);
			}
			MM(j, k) += sumInu;
		}
	}
}

for (int j = 2; j <= q; j++)
{
	for (int k = 1; k <= -1 + j; k++)
	{
		{
			scalar sumInu = 0;
			for (int nu = j; nu <= j; nu++)
			{
				sumInu += azero*vEG2G4c*b[nu] + azero*vEG2G4cAW*b[nu]*psiexp(-k + nu,1);
			}
			MM(j, k) += sumInu;
		}
	}
}

for (int j = 2; j <= q; j++)
{
	for (int k = 1; k <= -1 + j; k++)
	{
		{
			scalar sumInu = 0;
			for (int nu = k; nu <= k; nu++)
			{
				sumInu += azero*vEG4G2c*b[nu] + azero*vEG4G2cAD*b[nu]*psiexp(j - nu,1) + azero*vEG4G2cBD*b[nu]*psiexp(j - nu,2);
			}
			MM(j, k) += sumInu;
		}
	}
}

for (int j = 1; j <= -1 + q; j++)
{
	for (int k = j; k <= j; k++)
	{
		{
			scalar sumInu = 0;
			for (int nu = 1 + j; nu <= q; nu++)
			{
				sumInu += azero*vEG4G2c*b[nu] + azero*vEG4G2cAD*b[nu]*psiexp(-k + nu,1) + azero*vEG4G2cBD*b[nu]*psiexp(-k + nu,2);
			}
			MM(j, k) += sumInu;
		}
	}
}

for (int j = 1; j <= q; j++)
{
	for (int k = j; k <= j; k++)
	{
		{
			scalar sumInu = 0;
			for (int nu = 0; nu <= -1 + j; nu++)
			{
				sumInu += azero*vEG2G4c*b[nu] + azero*vEG2G4cAW*b[nu]*psiexp(k - nu,1);
			}
			MM(j, k) += sumInu;
		}
	}
}

for (int j = 1; j <= q; j++)
{
	for (int k = j; k <= j; k++)
	{
		{
			scalar sumInu = 0;
			for (int nu = j; nu <= j; nu++)
			{
				sumInu += azero*vEG6c*b[nu];
			}
			MM(j, k) += sumInu;
		}
	}
}

for (int j = 1; j <= -3 + q; j++)
{
	for (int k = 1 + j; k <= -2 + q; k++)
	{
		{
			scalar sumInuIkappa = 0;
			for (int nu = 1 + k; nu <= -1 + q; nu++)
			{
				for (int kappa = 1 + nu; kappa <= q; kappa++)
				{
					sumInuIkappa += vEG2G2G2G2c*b[kappa]*b[nu] + vEG2G2G2G2cAW*b[kappa]*b[nu]*psiexp(-j + k,1) + vEG2G2G2G2cAWADAH*b[kappa]*b[nu]*psiexp(-j + kappa,1) + vEG2G2G2G2cADAH*b[kappa]*b[nu]*psiexp(-k + kappa,1) + vEG2G2G2G2cAWBDBH*b[kappa]*b[nu]*psiexp(-j + k,1)*psiexp(-k + kappa,2) + vEG2G2G2G2cAH*b[kappa]*b[nu]*psiexp(kappa - nu,1) + vEG2G2G2G2cAWAH*b[kappa]*b[nu]*psiexp(-j + k + kappa - nu,1) + vEG2G2G2G2cAWAD*b[kappa]*b[nu]*psiexp(-j + nu,1) + vEG2G2G2G2cAWADBH*b[kappa]*b[nu]*psiexp(kappa - nu,2)*psiexp(-j + nu,1) + vEG2G2G2G2cAD*b[kappa]*b[nu]*psiexp(-k + nu,1) + vEG2G2G2G2cADBH*b[kappa]*b[nu]*psiexp(kappa - nu,2)*psiexp(-k + nu,1) + vEG2G2G2G2cAWBD*b[kappa]*b[nu]*psiexp(-j + k,1)*psiexp(-k + nu,2) + vEG2G2G2G2cAWBDCH*b[kappa]*b[nu]*psiexp(-j + k,1)*psiexp(kappa - nu,3)*psiexp(-k + nu,2) + vEG2G2G2G2cAWAHBD*b[kappa]*b[nu]*psiexp(-j + k + kappa - nu,1)*psiexp(-k + nu,2);
				}
			}
			MM(j, k) += sumInuIkappa;
		}
	}
}

for (int j = 1; j <= -3 + q; j++)
{
	for (int k = 1 + j; k <= -2 + q; k++)
	{
		{
			scalar sumInuIkappa = 0;
			for (int nu = 2 + k; nu <= q; nu++)
			{
				for (int kappa = 1 + k; kappa <= -1 + nu; kappa++)
				{
					sumInuIkappa += vEG2G2G2G2c*b[kappa]*b[nu] + vEG2G2G2G2cAW*b[kappa]*b[nu]*psiexp(-j + k,1) + vEG2G2G2G2cAWAD*b[kappa]*b[nu]*psiexp(-j + kappa,1) + vEG2G2G2G2cAD*b[kappa]*b[nu]*psiexp(-k + kappa,1) + vEG2G2G2G2cAWBD*b[kappa]*b[nu]*psiexp(-j + k,1)*psiexp(-k + kappa,2) + vEG2G2G2G2cAWADAH*b[kappa]*b[nu]*psiexp(-j + nu,1) + vEG2G2G2G2cADAH*b[kappa]*b[nu]*psiexp(-k + nu,1) + vEG2G2G2G2cAWBDBH*b[kappa]*b[nu]*psiexp(-j + k,1)*psiexp(-k + nu,2) + vEG2G2G2G2cAH*b[kappa]*b[nu]*psiexp(-kappa + nu,1) + vEG2G2G2G2cAWADBH*b[kappa]*b[nu]*psiexp(-j + kappa,1)*psiexp(-kappa + nu,2) + vEG2G2G2G2cADBH*b[kappa]*b[nu]*psiexp(-k + kappa,1)*psiexp(-kappa + nu,2) + vEG2G2G2G2cAWBDCH*b[kappa]*b[nu]*psiexp(-j + k,1)*psiexp(-k + kappa,2)*psiexp(-kappa + nu,3) + vEG2G2G2G2cAWAH*b[kappa]*b[nu]*psiexp(-j + k - kappa + nu,1) + vEG2G2G2G2cAWAHBD*b[kappa]*b[nu]*psiexp(-k + kappa,2)*psiexp(-j + k - kappa + nu,1);
				}
			}
			MM(j, k) += sumInuIkappa;
		}
	}
}

for (int j = 1; j <= -3 + q; j++)
{
	for (int k = 2 + j; k <= -1 + q; k++)
	{
		{
			scalar sumInuIkappa = 0;
			for (int nu = 1 + k; nu <= q; nu++)
			{
				for (int kappa = 1 + j; kappa <= -1 + k; kappa++)
				{
					sumInuIkappa += vEG2G2G2G2c*b[kappa]*b[nu] + vEG2G2G2G2cAWAD*b[kappa]*b[nu]*psiexp(-j + k,1) + vEG2G2G2G2cAD*b[kappa]*b[nu]*psiexp(k - kappa,1) + vEG2G2G2G2cAW*b[kappa]*b[nu]*psiexp(-j + kappa,1) + vEG2G2G2G2cAWBD*b[kappa]*b[nu]*psiexp(k - kappa,2)*psiexp(-j + kappa,1) + vEG2G2G2G2cAWADAH*b[kappa]*b[nu]*psiexp(-j + nu,1) + vEG2G2G2G2cAH*b[kappa]*b[nu]*psiexp(-k + nu,1) + vEG2G2G2G2cAWADBH*b[kappa]*b[nu]*psiexp(-j + k,1)*psiexp(-k + nu,2) + vEG2G2G2G2cADBH*b[kappa]*b[nu]*psiexp(k - kappa,1)*psiexp(-k + nu,2) + vEG2G2G2G2cAWBDCH*b[kappa]*b[nu]*psiexp(k - kappa,2)*psiexp(-j + kappa,1)*psiexp(-k + nu,3) + vEG2G2G2G2cADAH*b[kappa]*b[nu]*psiexp(-kappa + nu,1) + vEG2G2G2G2cAWBDBH*b[kappa]*b[nu]*psiexp(-j + kappa,1)*psiexp(-kappa + nu,2) + vEG2G2G2G2cAWAH*b[kappa]*b[nu]*psiexp(-j - k + kappa + nu,1) + vEG2G2G2G2cAWAHBD*b[kappa]*b[nu]*psiexp(k - kappa,2)*psiexp(-j - k + kappa + nu,1);
				}
			}
			MM(j, k) += sumInuIkappa;
		}
	}
}

for (int j = 1; j <= -2 + q; j++)
{
	for (int k = 1 + j; k <= -1 + q; k++)
	{
		{
			scalar sumInuIkappa = 0;
			for (int nu = 1 + k; nu <= q; nu++)
			{
				for (int kappa = 0; kappa <= -1 + j; kappa++)
				{
					sumInuIkappa += vEG2G2G2G2c*b[kappa]*b[nu] + vEG2G2G2G2cAD*b[kappa]*b[nu]*psiexp(-j + k,1) + vEG2G2G2G2cAW*b[kappa]*b[nu]*psiexp(j - kappa,1) + vEG2G2G2G2cAWBD*b[kappa]*b[nu]*psiexp(-j + k,2)*psiexp(j - kappa,1) + vEG2G2G2G2cAWAD*b[kappa]*b[nu]*psiexp(k - kappa,1) + vEG2G2G2G2cADAH*b[kappa]*b[nu]*psiexp(-j + nu,1) + vEG2G2G2G2cAWBDBH*b[kappa]*b[nu]*psiexp(j - kappa,1)*psiexp(-j + nu,2) + vEG2G2G2G2cAH*b[kappa]*b[nu]*psiexp(-k + nu,1) + vEG2G2G2G2cADBH*b[kappa]*b[nu]*psiexp(-j + k,1)*psiexp(-k + nu,2) + vEG2G2G2G2cAWADBH*b[kappa]*b[nu]*psiexp(k - kappa,1)*psiexp(-k + nu,2) + vEG2G2G2G2cAWBDCH*b[kappa]*b[nu]*psiexp(-j + k,2)*psiexp(j - kappa,1)*psiexp(-k + nu,3) + vEG2G2G2G2cAWADAH*b[kappa]*b[nu]*psiexp(-kappa + nu,1) + vEG2G2G2G2cAWAH*b[kappa]*b[nu]*psiexp(j - k - kappa + nu,1) + vEG2G2G2G2cAWAHBD*b[kappa]*b[nu]*psiexp(-j + k,2)*psiexp(j - k - kappa + nu,1);
				}
			}
			MM(j, k) += sumInuIkappa;
		}
	}
}

for (int j = 1; j <= -2 + q; j++)
{
	for (int k = 1 + j; k <= -1 + q; k++)
	{
		{
			scalar sumInuIkappa = 0;
			for (int nu = 1 + k; nu <= q; nu++)
			{
				for (int kappa = nu; kappa <= nu; kappa++)
				{
					sumInuIkappa += vEG2G2G4c*Power2(b[nu]) + vEG2G2G4cAW*Power2(b[nu])*psiexp(-j + k,1) + vEG2G2G4cAWAD*Power2(b[nu])*psiexp(-j + nu,1) + vEG2G2G4cAD*Power2(b[nu])*psiexp(-k + nu,1) + vEG2G2G4cAWBD*Power2(b[nu])*psiexp(-j + k,1)*psiexp(-k + nu,2);
				}
			}
			MM(j, k) += sumInuIkappa;
		}
	}
}

for (int j = 1; j <= -2 + q; j++)
{
	for (int k = 1 + j; k <= -1 + q; k++)
	{
		{
			scalar sumInuIkappa = 0;
			for (int nu = 1 + k; nu <= q; nu++)
			{
				for (int kappa = k; kappa <= k; kappa++)
				{
					sumInuIkappa += vEG2G4G2c*b[kappa]*b[nu] + vEG2G4G2cAW*b[kappa]*b[nu]*psiexp(-j + kappa,1) + vEG2G4G2cAWAH*b[kappa]*b[nu]*psiexp(-j + nu,1) + vEG2G4G2cAH*b[kappa]*b[nu]*psiexp(-kappa + nu,1) + vEG2G4G2cBH*b[kappa]*b[nu]*psiexp(-kappa + nu,2) + vEG2G4G2cAWBH*b[kappa]*b[nu]*psiexp(-j + kappa,1)*psiexp(-kappa + nu,2) + vEG2G4G2cAWCH*b[kappa]*b[nu]*psiexp(-j + kappa,1)*psiexp(-kappa + nu,3);
				}
			}
			MM(j, k) += sumInuIkappa;
		}
	}
}

for (int j = 1; j <= -2 + q; j++)
{
	for (int k = 1 + j; k <= -1 + q; k++)
	{
		{
			scalar sumInuIkappa = 0;
			for (int nu = 1 + k; nu <= q; nu++)
			{
				for (int kappa = j; kappa <= j; kappa++)
				{
					sumInuIkappa += vEG4G2G2c*b[kappa]*b[nu] + vEG4G2G2cAD*b[kappa]*b[nu]*psiexp(k - kappa,1) + vEG4G2G2cBD*b[kappa]*b[nu]*psiexp(k - kappa,2) + vEG4G2G2cAH*b[kappa]*b[nu]*psiexp(-k + nu,1) + vEG4G2G2cAHBD*b[kappa]*b[nu]*psiexp(k - kappa,2)*psiexp(-k + nu,1) + vEG4G2G2cADBH*b[kappa]*b[nu]*psiexp(k - kappa,1)*psiexp(-k + nu,2) + vEG4G2G2cBDCH*b[kappa]*b[nu]*psiexp(k - kappa,2)*psiexp(-k + nu,3) + vEG4G2G2cADAH*b[kappa]*b[nu]*psiexp(-kappa + nu,1) + vEG4G2G2cBDBH*b[kappa]*b[nu]*psiexp(-kappa + nu,2);
				}
			}
			MM(j, k) += sumInuIkappa;
		}
	}
}

for (int j = 1; j <= -3 + q; j++)
{
	for (int k = 2 + j; k <= -1 + q; k++)
	{
		{
			scalar sumInuIkappa = 0;
			for (int nu = 1 + j; nu <= -1 + k; nu++)
			{
				for (int kappa = 1 + k; kappa <= q; kappa++)
				{
					sumInuIkappa += vEG2G2G2G2c*b[kappa]*b[nu] + vEG2G2G2G2cAWAD*b[kappa]*b[nu]*psiexp(-j + k,1) + vEG2G2G2G2cAWADAH*b[kappa]*b[nu]*psiexp(-j + kappa,1) + vEG2G2G2G2cAH*b[kappa]*b[nu]*psiexp(-k + kappa,1) + vEG2G2G2G2cAWADBH*b[kappa]*b[nu]*psiexp(-j + k,1)*psiexp(-k + kappa,2) + vEG2G2G2G2cAD*b[kappa]*b[nu]*psiexp(k - nu,1) + vEG2G2G2G2cADBH*b[kappa]*b[nu]*psiexp(-k + kappa,2)*psiexp(k - nu,1) + vEG2G2G2G2cADAH*b[kappa]*b[nu]*psiexp(kappa - nu,1) + vEG2G2G2G2cAW*b[kappa]*b[nu]*psiexp(-j + nu,1) + vEG2G2G2G2cAWBD*b[kappa]*b[nu]*psiexp(k - nu,2)*psiexp(-j + nu,1) + vEG2G2G2G2cAWBDCH*b[kappa]*b[nu]*psiexp(-k + kappa,3)*psiexp(k - nu,2)*psiexp(-j + nu,1) + vEG2G2G2G2cAWBDBH*b[kappa]*b[nu]*psiexp(kappa - nu,2)*psiexp(-j + nu,1) + vEG2G2G2G2cAWAH*b[kappa]*b[nu]*psiexp(-j - k + kappa + nu,1) + vEG2G2G2G2cAWAHBD*b[kappa]*b[nu]*psiexp(k - nu,2)*psiexp(-j - k + kappa + nu,1);
				}
			}
			MM(j, k) += sumInuIkappa;
		}
	}
}

for (int j = 1; j <= -3 + q; j++)
{
	for (int k = 3 + j; k <= q; k++)
	{
		{
			scalar sumInuIkappa = 0;
			for (int nu = 1 + j; nu <= -2 + k; nu++)
			{
				for (int kappa = 1 + nu; kappa <= -1 + k; kappa++)
				{
					sumInuIkappa += vEG2G2G2G2c*b[kappa]*b[nu] + vEG2G2G2G2cAWADAH*b[kappa]*b[nu]*psiexp(-j + k,1) + vEG2G2G2G2cAH*b[kappa]*b[nu]*psiexp(k - kappa,1) + vEG2G2G2G2cAWAD*b[kappa]*b[nu]*psiexp(-j + kappa,1) + vEG2G2G2G2cAWADBH*b[kappa]*b[nu]*psiexp(k - kappa,2)*psiexp(-j + kappa,1) + vEG2G2G2G2cADAH*b[kappa]*b[nu]*psiexp(k - nu,1) + vEG2G2G2G2cAD*b[kappa]*b[nu]*psiexp(kappa - nu,1) + vEG2G2G2G2cADBH*b[kappa]*b[nu]*psiexp(k - kappa,2)*psiexp(kappa - nu,1) + vEG2G2G2G2cAW*b[kappa]*b[nu]*psiexp(-j + nu,1) + vEG2G2G2G2cAWBDBH*b[kappa]*b[nu]*psiexp(k - nu,2)*psiexp(-j + nu,1) + vEG2G2G2G2cAWBD*b[kappa]*b[nu]*psiexp(kappa - nu,2)*psiexp(-j + nu,1) + vEG2G2G2G2cAWBDCH*b[kappa]*b[nu]*psiexp(k - kappa,3)*psiexp(kappa - nu,2)*psiexp(-j + nu,1) + vEG2G2G2G2cAWAH*b[kappa]*b[nu]*psiexp(-j + k - kappa + nu,1) + vEG2G2G2G2cAWAHBD*b[kappa]*b[nu]*psiexp(kappa - nu,2)*psiexp(-j + k - kappa + nu,1);
				}
			}
			MM(j, k) += sumInuIkappa;
		}
	}
}

for (int j = 1; j <= -3 + q; j++)
{
	for (int k = 3 + j; k <= q; k++)
	{
		{
			scalar sumInuIkappa = 0;
			for (int nu = 2 + j; nu <= -1 + k; nu++)
			{
				for (int kappa = 1 + j; kappa <= -1 + nu; kappa++)
				{
					sumInuIkappa += vEG2G2G2G2c*b[kappa]*b[nu] + vEG2G2G2G2cAWADAH*b[kappa]*b[nu]*psiexp(-j + k,1) + vEG2G2G2G2cADAH*b[kappa]*b[nu]*psiexp(k - kappa,1) + vEG2G2G2G2cAW*b[kappa]*b[nu]*psiexp(-j + kappa,1) + vEG2G2G2G2cAWBDBH*b[kappa]*b[nu]*psiexp(k - kappa,2)*psiexp(-j + kappa,1) + vEG2G2G2G2cAH*b[kappa]*b[nu]*psiexp(k - nu,1) + vEG2G2G2G2cAWAH*b[kappa]*b[nu]*psiexp(-j + k + kappa - nu,1) + vEG2G2G2G2cAWAD*b[kappa]*b[nu]*psiexp(-j + nu,1) + vEG2G2G2G2cAWADBH*b[kappa]*b[nu]*psiexp(k - nu,2)*psiexp(-j + nu,1) + vEG2G2G2G2cAD*b[kappa]*b[nu]*psiexp(-kappa + nu,1) + vEG2G2G2G2cADBH*b[kappa]*b[nu]*psiexp(k - nu,2)*psiexp(-kappa + nu,1) + vEG2G2G2G2cAWBD*b[kappa]*b[nu]*psiexp(-j + kappa,1)*psiexp(-kappa + nu,2) + vEG2G2G2G2cAWBDCH*b[kappa]*b[nu]*psiexp(-j + kappa,1)*psiexp(k - nu,3)*psiexp(-kappa + nu,2) + vEG2G2G2G2cAWAHBD*b[kappa]*b[nu]*psiexp(-j + k + kappa - nu,1)*psiexp(-kappa + nu,2);
				}
			}
			MM(j, k) += sumInuIkappa;
		}
	}
}

for (int j = 1; j <= -2 + q; j++)
{
	for (int k = 2 + j; k <= q; k++)
	{
		{
			scalar sumInuIkappa = 0;
			for (int nu = 1 + j; nu <= -1 + k; nu++)
			{
				for (int kappa = 0; kappa <= -1 + j; kappa++)
				{
					sumInuIkappa += vEG2G2G2G2c*b[kappa]*b[nu] + vEG2G2G2G2cADAH*b[kappa]*b[nu]*psiexp(-j + k,1) + vEG2G2G2G2cAW*b[kappa]*b[nu]*psiexp(j - kappa,1) + vEG2G2G2G2cAWBDBH*b[kappa]*b[nu]*psiexp(-j + k,2)*psiexp(j - kappa,1) + vEG2G2G2G2cAWADAH*b[kappa]*b[nu]*psiexp(k - kappa,1) + vEG2G2G2G2cAH*b[kappa]*b[nu]*psiexp(k - nu,1) + vEG2G2G2G2cAWAH*b[kappa]*b[nu]*psiexp(j + k - kappa - nu,1) + vEG2G2G2G2cAD*b[kappa]*b[nu]*psiexp(-j + nu,1) + vEG2G2G2G2cADBH*b[kappa]*b[nu]*psiexp(k - nu,2)*psiexp(-j + nu,1) + vEG2G2G2G2cAWBD*b[kappa]*b[nu]*psiexp(j - kappa,1)*psiexp(-j + nu,2) + vEG2G2G2G2cAWBDCH*b[kappa]*b[nu]*psiexp(j - kappa,1)*psiexp(k - nu,3)*psiexp(-j + nu,2) + vEG2G2G2G2cAWAHBD*b[kappa]*b[nu]*psiexp(j + k - kappa - nu,1)*psiexp(-j + nu,2) + vEG2G2G2G2cAWAD*b[kappa]*b[nu]*psiexp(-kappa + nu,1) + vEG2G2G2G2cAWADBH*b[kappa]*b[nu]*psiexp(k - nu,2)*psiexp(-kappa + nu,1);
				}
			}
			MM(j, k) += sumInuIkappa;
		}
	}
}

for (int j = 1; j <= -2 + q; j++)
{
	for (int k = 2 + j; k <= q; k++)
	{
		{
			scalar sumInuIkappa = 0;
			for (int nu = 1 + j; nu <= -1 + k; nu++)
			{
				for (int kappa = k; kappa <= k; kappa++)
				{
					sumInuIkappa += vEG2G2G4c*b[kappa]*b[nu] + vEG2G2G4cAWAD*b[kappa]*b[nu]*psiexp(-j + kappa,1) + vEG2G2G4cAD*b[kappa]*b[nu]*psiexp(kappa - nu,1) + vEG2G2G4cAW*b[kappa]*b[nu]*psiexp(-j + nu,1) + vEG2G2G4cAWBD*b[kappa]*b[nu]*psiexp(kappa - nu,2)*psiexp(-j + nu,1);
				}
			}
			MM(j, k) += sumInuIkappa;
		}
	}
}

for (int j = 1; j <= -2 + q; j++)
{
	for (int k = 2 + j; k <= q; k++)
	{
		{
			scalar sumInuIkappa = 0;
			for (int nu = 1 + j; nu <= -1 + k; nu++)
			{
				for (int kappa = nu; kappa <= nu; kappa++)
				{
					sumInuIkappa += vEG2G4G2c*Power2(b[nu]) + vEG2G4G2cAWAH*Power2(b[nu])*psiexp(-j + k,1) + vEG2G4G2cAH*Power2(b[nu])*psiexp(k - nu,1) + vEG2G4G2cBH*Power2(b[nu])*psiexp(k - nu,2) + vEG2G4G2cAW*Power2(b[nu])*psiexp(-j + nu,1) + vEG2G4G2cAWBH*Power2(b[nu])*psiexp(k - nu,2)*psiexp(-j + nu,1) + vEG2G4G2cAWCH*Power2(b[nu])*psiexp(k - nu,3)*psiexp(-j + nu,1);
				}
			}
			MM(j, k) += sumInuIkappa;
		}
	}
}

for (int j = 1; j <= -2 + q; j++)
{
	for (int k = 2 + j; k <= q; k++)
	{
		{
			scalar sumInuIkappa = 0;
			for (int nu = 1 + j; nu <= -1 + k; nu++)
			{
				for (int kappa = j; kappa <= j; kappa++)
				{
					sumInuIkappa += vEG4G2G2c*b[kappa]*b[nu] + vEG4G2G2cADAH*b[kappa]*b[nu]*psiexp(k - kappa,1) + vEG4G2G2cBDBH*b[kappa]*b[nu]*psiexp(k - kappa,2) + vEG4G2G2cAH*b[kappa]*b[nu]*psiexp(k - nu,1) + vEG4G2G2cAD*b[kappa]*b[nu]*psiexp(-kappa + nu,1) + vEG4G2G2cADBH*b[kappa]*b[nu]*psiexp(k - nu,2)*psiexp(-kappa + nu,1) + vEG4G2G2cBD*b[kappa]*b[nu]*psiexp(-kappa + nu,2) + vEG4G2G2cAHBD*b[kappa]*b[nu]*psiexp(k - nu,1)*psiexp(-kappa + nu,2) + vEG4G2G2cBDCH*b[kappa]*b[nu]*psiexp(k - nu,3)*psiexp(-kappa + nu,2);
				}
			}
			MM(j, k) += sumInuIkappa;
		}
	}
}

for (int j = 1; j <= -2 + q; j++)
{
	for (int k = 1 + j; k <= -1 + q; k++)
	{
		{
			scalar sumInuIkappa = 0;
			for (int nu = 0; nu <= -1 + j; nu++)
			{
				for (int kappa = 1 + k; kappa <= q; kappa++)
				{
					sumInuIkappa += vEG2G2G2G2c*b[kappa]*b[nu] + vEG2G2G2G2cAD*b[kappa]*b[nu]*psiexp(-j + k,1) + vEG2G2G2G2cADAH*b[kappa]*b[nu]*psiexp(-j + kappa,1) + vEG2G2G2G2cAH*b[kappa]*b[nu]*psiexp(-k + kappa,1) + vEG2G2G2G2cADBH*b[kappa]*b[nu]*psiexp(-j + k,1)*psiexp(-k + kappa,2) + vEG2G2G2G2cAW*b[kappa]*b[nu]*psiexp(j - nu,1) + vEG2G2G2G2cAWBD*b[kappa]*b[nu]*psiexp(-j + k,2)*psiexp(j - nu,1) + vEG2G2G2G2cAWBDBH*b[kappa]*b[nu]*psiexp(-j + kappa,2)*psiexp(j - nu,1) + vEG2G2G2G2cAWBDCH*b[kappa]*b[nu]*psiexp(-j + k,2)*psiexp(-k + kappa,3)*psiexp(j - nu,1) + vEG2G2G2G2cAWAD*b[kappa]*b[nu]*psiexp(k - nu,1) + vEG2G2G2G2cAWADBH*b[kappa]*b[nu]*psiexp(-k + kappa,2)*psiexp(k - nu,1) + vEG2G2G2G2cAWADAH*b[kappa]*b[nu]*psiexp(kappa - nu,1) + vEG2G2G2G2cAWAH*b[kappa]*b[nu]*psiexp(j - k + kappa - nu,1) + vEG2G2G2G2cAWAHBD*b[kappa]*b[nu]*psiexp(-j + k,2)*psiexp(j - k + kappa - nu,1);
				}
			}
			MM(j, k) += sumInuIkappa;
		}
	}
}

for (int j = 1; j <= -2 + q; j++)
{
	for (int k = 2 + j; k <= q; k++)
	{
		{
			scalar sumInuIkappa = 0;
			for (int nu = 0; nu <= -1 + j; nu++)
			{
				for (int kappa = 1 + j; kappa <= -1 + k; kappa++)
				{
					sumInuIkappa += vEG2G2G2G2c*b[kappa]*b[nu] + vEG2G2G2G2cADAH*b[kappa]*b[nu]*psiexp(-j + k,1) + vEG2G2G2G2cAH*b[kappa]*b[nu]*psiexp(k - kappa,1) + vEG2G2G2G2cAD*b[kappa]*b[nu]*psiexp(-j + kappa,1) + vEG2G2G2G2cADBH*b[kappa]*b[nu]*psiexp(k - kappa,2)*psiexp(-j + kappa,1) + vEG2G2G2G2cAW*b[kappa]*b[nu]*psiexp(j - nu,1) + vEG2G2G2G2cAWBDBH*b[kappa]*b[nu]*psiexp(-j + k,2)*psiexp(j - nu,1) + vEG2G2G2G2cAWBD*b[kappa]*b[nu]*psiexp(-j + kappa,2)*psiexp(j - nu,1) + vEG2G2G2G2cAWBDCH*b[kappa]*b[nu]*psiexp(k - kappa,3)*psiexp(-j + kappa,2)*psiexp(j - nu,1) + vEG2G2G2G2cAWADAH*b[kappa]*b[nu]*psiexp(k - nu,1) + vEG2G2G2G2cAWAH*b[kappa]*b[nu]*psiexp(j + k - kappa - nu,1) + vEG2G2G2G2cAWAHBD*b[kappa]*b[nu]*psiexp(-j + kappa,2)*psiexp(j + k - kappa - nu,1) + vEG2G2G2G2cAWAD*b[kappa]*b[nu]*psiexp(kappa - nu,1) + vEG2G2G2G2cAWADBH*b[kappa]*b[nu]*psiexp(k - kappa,2)*psiexp(kappa - nu,1);
				}
			}
			MM(j, k) += sumInuIkappa;
		}
	}
}

for (int j = 2; j <= -1 + q; j++)
{
	for (int k = 1 + j; k <= q; k++)
	{
		{
			scalar sumInuIkappa = 0;
			for (int nu = 0; nu <= -2 + j; nu++)
			{
				for (int kappa = 1 + nu; kappa <= -1 + j; kappa++)
				{
					sumInuIkappa += vEG2G2G2G2c*b[kappa]*b[nu] + vEG2G2G2G2cAH*b[kappa]*b[nu]*psiexp(-j + k,1) + vEG2G2G2G2cAD*b[kappa]*b[nu]*psiexp(j - kappa,1) + vEG2G2G2G2cADBH*b[kappa]*b[nu]*psiexp(-j + k,2)*psiexp(j - kappa,1) + vEG2G2G2G2cADAH*b[kappa]*b[nu]*psiexp(k - kappa,1) + vEG2G2G2G2cAWAD*b[kappa]*b[nu]*psiexp(j - nu,1) + vEG2G2G2G2cAWADBH*b[kappa]*b[nu]*psiexp(-j + k,2)*psiexp(j - nu,1) + vEG2G2G2G2cAWADAH*b[kappa]*b[nu]*psiexp(k - nu,1) + vEG2G2G2G2cAW*b[kappa]*b[nu]*psiexp(kappa - nu,1) + vEG2G2G2G2cAWBD*b[kappa]*b[nu]*psiexp(j - kappa,2)*psiexp(kappa - nu,1) + vEG2G2G2G2cAWBDCH*b[kappa]*b[nu]*psiexp(-j + k,3)*psiexp(j - kappa,2)*psiexp(kappa - nu,1) + vEG2G2G2G2cAWBDBH*b[kappa]*b[nu]*psiexp(k - kappa,2)*psiexp(kappa - nu,1) + vEG2G2G2G2cAWAH*b[kappa]*b[nu]*psiexp(-j + k + kappa - nu,1) + vEG2G2G2G2cAWAHBD*b[kappa]*b[nu]*psiexp(j - kappa,2)*psiexp(-j + k + kappa - nu,1);
				}
			}
			MM(j, k) += sumInuIkappa;
		}
	}
}

for (int j = 2; j <= -1 + q; j++)
{
	for (int k = 1 + j; k <= q; k++)
	{
		{
			scalar sumInuIkappa = 0;
			for (int nu = 1; nu <= -1 + j; nu++)
			{
				for (int kappa = 0; kappa <= -1 + nu; kappa++)
				{
					sumInuIkappa += vEG2G2G2G2c*b[kappa]*b[nu] + vEG2G2G2G2cAH*b[kappa]*b[nu]*psiexp(-j + k,1) + vEG2G2G2G2cAWAD*b[kappa]*b[nu]*psiexp(j - kappa,1) + vEG2G2G2G2cAWADBH*b[kappa]*b[nu]*psiexp(-j + k,2)*psiexp(j - kappa,1) + vEG2G2G2G2cAWADAH*b[kappa]*b[nu]*psiexp(k - kappa,1) + vEG2G2G2G2cAD*b[kappa]*b[nu]*psiexp(j - nu,1) + vEG2G2G2G2cADBH*b[kappa]*b[nu]*psiexp(-j + k,2)*psiexp(j - nu,1) + vEG2G2G2G2cADAH*b[kappa]*b[nu]*psiexp(k - nu,1) + vEG2G2G2G2cAW*b[kappa]*b[nu]*psiexp(-kappa + nu,1) + vEG2G2G2G2cAWBD*b[kappa]*b[nu]*psiexp(j - nu,2)*psiexp(-kappa + nu,1) + vEG2G2G2G2cAWBDCH*b[kappa]*b[nu]*psiexp(-j + k,3)*psiexp(j - nu,2)*psiexp(-kappa + nu,1) + vEG2G2G2G2cAWBDBH*b[kappa]*b[nu]*psiexp(k - nu,2)*psiexp(-kappa + nu,1) + vEG2G2G2G2cAWAH*b[kappa]*b[nu]*psiexp(-j + k - kappa + nu,1) + vEG2G2G2G2cAWAHBD*b[kappa]*b[nu]*psiexp(j - nu,2)*psiexp(-j + k - kappa + nu,1);
				}
			}
			MM(j, k) += sumInuIkappa;
		}
	}
}

for (int j = 1; j <= -1 + q; j++)
{
	for (int k = 1 + j; k <= q; k++)
	{
		{
			scalar sumInuIkappa = 0;
			for (int nu = 0; nu <= -1 + j; nu++)
			{
				for (int kappa = k; kappa <= k; kappa++)
				{
					sumInuIkappa += vEG2G2G4c*b[kappa]*b[nu] + vEG2G2G4cAD*b[kappa]*b[nu]*psiexp(-j + kappa,1) + vEG2G2G4cAW*b[kappa]*b[nu]*psiexp(j - nu,1) + vEG2G2G4cAWBD*b[kappa]*b[nu]*psiexp(-j + kappa,2)*psiexp(j - nu,1) + vEG2G2G4cAWAD*b[kappa]*b[nu]*psiexp(kappa - nu,1);
				}
			}
			MM(j, k) += sumInuIkappa;
		}
	}
}

for (int j = 1; j <= -1 + q; j++)
{
	for (int k = 1 + j; k <= q; k++)
	{
		{
			scalar sumInuIkappa = 0;
			for (int nu = 0; nu <= -1 + j; nu++)
			{
				for (int kappa = j; kappa <= j; kappa++)
				{
					sumInuIkappa += vEG2G4G2c*b[kappa]*b[nu] + vEG2G4G2cAH*b[kappa]*b[nu]*psiexp(k - kappa,1) + vEG2G4G2cBH*b[kappa]*b[nu]*psiexp(k - kappa,2) + vEG2G4G2cAWAH*b[kappa]*b[nu]*psiexp(k - nu,1) + vEG2G4G2cAW*b[kappa]*b[nu]*psiexp(kappa - nu,1) + vEG2G4G2cAWBH*b[kappa]*b[nu]*psiexp(k - kappa,2)*psiexp(kappa - nu,1) + vEG2G4G2cAWCH*b[kappa]*b[nu]*psiexp(k - kappa,3)*psiexp(kappa - nu,1);
				}
			}
			MM(j, k) += sumInuIkappa;
		}
	}
}

for (int j = 1; j <= -1 + q; j++)
{
	for (int k = 1 + j; k <= q; k++)
	{
		{
			scalar sumInuIkappa = 0;
			for (int nu = 0; nu <= -1 + j; nu++)
			{
				for (int kappa = nu; kappa <= nu; kappa++)
				{
					sumInuIkappa += vEG4G2G2c*Power2(b[nu]) + vEG4G2G2cAH*Power2(b[nu])*psiexp(-j + k,1) + vEG4G2G2cAD*Power2(b[nu])*psiexp(j - nu,1) + vEG4G2G2cADBH*Power2(b[nu])*psiexp(-j + k,2)*psiexp(j - nu,1) + vEG4G2G2cBD*Power2(b[nu])*psiexp(j - nu,2) + vEG4G2G2cAHBD*Power2(b[nu])*psiexp(-j + k,1)*psiexp(j - nu,2) + vEG4G2G2cBDCH*Power2(b[nu])*psiexp(-j + k,3)*psiexp(j - nu,2) + vEG4G2G2cADAH*Power2(b[nu])*psiexp(k - nu,1) + vEG4G2G2cBDBH*Power2(b[nu])*psiexp(k - nu,2);
				}
			}
			MM(j, k) += sumInuIkappa;
		}
	}
}

for (int j = 1; j <= -2 + q; j++)
{
	for (int k = 1 + j; k <= -1 + q; k++)
	{
		{
			scalar sumInuIkappa = 0;
			for (int nu = k; nu <= k; nu++)
			{
				for (int kappa = 1 + k; kappa <= q; kappa++)
				{
					sumInuIkappa += vEG2G4G2c*b[kappa]*b[nu] + vEG2G4G2cAWAH*b[kappa]*b[nu]*psiexp(-j + kappa,1) + vEG2G4G2cAH*b[kappa]*b[nu]*psiexp(kappa - nu,1) + vEG2G4G2cBH*b[kappa]*b[nu]*psiexp(kappa - nu,2) + vEG2G4G2cAW*b[kappa]*b[nu]*psiexp(-j + nu,1) + vEG2G4G2cAWBH*b[kappa]*b[nu]*psiexp(kappa - nu,2)*psiexp(-j + nu,1) + vEG2G4G2cAWCH*b[kappa]*b[nu]*psiexp(kappa - nu,3)*psiexp(-j + nu,1);
				}
			}
			MM(j, k) += sumInuIkappa;
		}
	}
}

for (int j = 1; j <= -2 + q; j++)
{
	for (int k = 2 + j; k <= q; k++)
	{
		{
			scalar sumInuIkappa = 0;
			for (int nu = k; nu <= k; nu++)
			{
				for (int kappa = 1 + j; kappa <= -1 + k; kappa++)
				{
					sumInuIkappa += vEG2G2G4c*b[kappa]*b[nu] + vEG2G2G4cAW*b[kappa]*b[nu]*psiexp(-j + kappa,1) + vEG2G2G4cAWAD*b[kappa]*b[nu]*psiexp(-j + nu,1) + vEG2G2G4cAD*b[kappa]*b[nu]*psiexp(-kappa + nu,1) + vEG2G2G4cAWBD*b[kappa]*b[nu]*psiexp(-j + kappa,1)*psiexp(-kappa + nu,2);
				}
			}
			MM(j, k) += sumInuIkappa;
		}
	}
}

for (int j = 1; j <= -1 + q; j++)
{
	for (int k = 1 + j; k <= q; k++)
	{
		{
			scalar sumInuIkappa = 0;
			for (int nu = k; nu <= k; nu++)
			{
				for (int kappa = 0; kappa <= -1 + j; kappa++)
				{
					sumInuIkappa += vEG2G2G4c*b[kappa]*b[nu] + vEG2G2G4cAW*b[kappa]*b[nu]*psiexp(j - kappa,1) + vEG2G2G4cAD*b[kappa]*b[nu]*psiexp(-j + nu,1) + vEG2G2G4cAWBD*b[kappa]*b[nu]*psiexp(j - kappa,1)*psiexp(-j + nu,2) + vEG2G2G4cAWAD*b[kappa]*b[nu]*psiexp(-kappa + nu,1);
				}
			}
			MM(j, k) += sumInuIkappa;
		}
	}
}

for (int j = 1; j <= -1 + q; j++)
{
	for (int k = 1 + j; k <= q; k++)
	{
		{
			scalar sumInuIkappa = 0;
			for (int nu = k; nu <= k; nu++)
			{
				for (int kappa = k; kappa <= k; kappa++)
				{
					sumInuIkappa += vEG2G6c*Power2(b[nu]) + vEG2G6cAW*Power2(b[nu])*psiexp(-j + nu,1);
				}
			}
			MM(j, k) += sumInuIkappa;
		}
	}
}

for (int j = 1; j <= -1 + q; j++)
{
	for (int k = 1 + j; k <= q; k++)
	{
		{
			scalar sumInuIkappa = 0;
			for (int nu = k; nu <= k; nu++)
			{
				for (int kappa = j; kappa <= j; kappa++)
				{
					sumInuIkappa += vEG4G4c*b[kappa]*b[nu] + vEG4G4cAD*b[kappa]*b[nu]*psiexp(-kappa + nu,1) + vEG4G4cBD*b[kappa]*b[nu]*psiexp(-kappa + nu,2);
				}
			}
			MM(j, k) += sumInuIkappa;
		}
	}
}

for (int j = 1; j <= -2 + q; j++)
{
	for (int k = 1 + j; k <= -1 + q; k++)
	{
		{
			scalar sumInuIkappa = 0;
			for (int nu = j; nu <= j; nu++)
			{
				for (int kappa = 1 + k; kappa <= q; kappa++)
				{
					sumInuIkappa += vEG4G2G2c*b[kappa]*b[nu] + vEG4G2G2cAH*b[kappa]*b[nu]*psiexp(-k + kappa,1) + vEG4G2G2cAD*b[kappa]*b[nu]*psiexp(k - nu,1) + vEG4G2G2cADBH*b[kappa]*b[nu]*psiexp(-k + kappa,2)*psiexp(k - nu,1) + vEG4G2G2cBD*b[kappa]*b[nu]*psiexp(k - nu,2) + vEG4G2G2cAHBD*b[kappa]*b[nu]*psiexp(-k + kappa,1)*psiexp(k - nu,2) + vEG4G2G2cBDCH*b[kappa]*b[nu]*psiexp(-k + kappa,3)*psiexp(k - nu,2) + vEG4G2G2cADAH*b[kappa]*b[nu]*psiexp(kappa - nu,1) + vEG4G2G2cBDBH*b[kappa]*b[nu]*psiexp(kappa - nu,2);
				}
			}
			MM(j, k) += sumInuIkappa;
		}
	}
}

for (int j = 1; j <= -2 + q; j++)
{
	for (int k = 2 + j; k <= q; k++)
	{
		{
			scalar sumInuIkappa = 0;
			for (int nu = j; nu <= j; nu++)
			{
				for (int kappa = 1 + j; kappa <= -1 + k; kappa++)
				{
					sumInuIkappa += vEG4G2G2c*b[kappa]*b[nu] + vEG4G2G2cAH*b[kappa]*b[nu]*psiexp(k - kappa,1) + vEG4G2G2cADAH*b[kappa]*b[nu]*psiexp(k - nu,1) + vEG4G2G2cBDBH*b[kappa]*b[nu]*psiexp(k - nu,2) + vEG4G2G2cAD*b[kappa]*b[nu]*psiexp(kappa - nu,1) + vEG4G2G2cADBH*b[kappa]*b[nu]*psiexp(k - kappa,2)*psiexp(kappa - nu,1) + vEG4G2G2cBD*b[kappa]*b[nu]*psiexp(kappa - nu,2) + vEG4G2G2cAHBD*b[kappa]*b[nu]*psiexp(k - kappa,1)*psiexp(kappa - nu,2) + vEG4G2G2cBDCH*b[kappa]*b[nu]*psiexp(k - kappa,3)*psiexp(kappa - nu,2);
				}
			}
			MM(j, k) += sumInuIkappa;
		}
	}
}

for (int j = 1; j <= -1 + q; j++)
{
	for (int k = 1 + j; k <= q; k++)
	{
		{
			scalar sumInuIkappa = 0;
			for (int nu = j; nu <= j; nu++)
			{
				for (int kappa = 0; kappa <= -1 + j; kappa++)
				{
					sumInuIkappa += vEG2G4G2c*b[kappa]*b[nu] + vEG2G4G2cAWAH*b[kappa]*b[nu]*psiexp(k - kappa,1) + vEG2G4G2cAH*b[kappa]*b[nu]*psiexp(k - nu,1) + vEG2G4G2cBH*b[kappa]*b[nu]*psiexp(k - nu,2) + vEG2G4G2cAW*b[kappa]*b[nu]*psiexp(-kappa + nu,1) + vEG2G4G2cAWBH*b[kappa]*b[nu]*psiexp(k - nu,2)*psiexp(-kappa + nu,1) + vEG2G4G2cAWCH*b[kappa]*b[nu]*psiexp(k - nu,3)*psiexp(-kappa + nu,1);
				}
			}
			MM(j, k) += sumInuIkappa;
		}
	}
}

for (int j = 1; j <= -1 + q; j++)
{
	for (int k = 1 + j; k <= q; k++)
	{
		{
			scalar sumInuIkappa = 0;
			for (int nu = j; nu <= j; nu++)
			{
				for (int kappa = k; kappa <= k; kappa++)
				{
					sumInuIkappa += vEG4G4c*b[kappa]*b[nu] + vEG4G4cAD*b[kappa]*b[nu]*psiexp(kappa - nu,1) + vEG4G4cBD*b[kappa]*b[nu]*psiexp(kappa - nu,2);
				}
			}
			MM(j, k) += sumInuIkappa;
		}
	}
}

for (int j = 1; j <= -1 + q; j++)
{
	for (int k = 1 + j; k <= q; k++)
	{
		{
			scalar sumInuIkappa = 0;
			for (int nu = j; nu <= j; nu++)
			{
				for (int kappa = j; kappa <= j; kappa++)
				{
					sumInuIkappa += vEG6G2c*Power2(b[nu]) + vEG6G2cAH*Power2(b[nu])*psiexp(k - nu,1) + vEG6G2cBH*Power2(b[nu])*psiexp(k - nu,2) + vEG6G2cCH*Power2(b[nu])*psiexp(k - nu,3);
				}
			}
			MM(j, k) += sumInuIkappa;
		}
	}
}

for (int j = 2; j <= -2 + q; j++)
{
	for (int k = 1; k <= -1 + j; k++)
	{
		{
			scalar sumInuIkappa = 0;
			for (int nu = 1 + j; nu <= -1 + q; nu++)
			{
				for (int kappa = 1 + nu; kappa <= q; kappa++)
				{
					sumInuIkappa += vEG2G2G2G2c*b[kappa]*b[nu] + vEG2G2G2G2cAW*b[kappa]*b[nu]*psiexp(j - k,1) + vEG2G2G2G2cADAH*b[kappa]*b[nu]*psiexp(-j + kappa,1) + vEG2G2G2G2cAWBDBH*b[kappa]*b[nu]*psiexp(j - k,1)*psiexp(-j + kappa,2) + vEG2G2G2G2cAWADAH*b[kappa]*b[nu]*psiexp(-k + kappa,1) + vEG2G2G2G2cAH*b[kappa]*b[nu]*psiexp(kappa - nu,1) + vEG2G2G2G2cAWAH*b[kappa]*b[nu]*psiexp(j - k + kappa - nu,1) + vEG2G2G2G2cAD*b[kappa]*b[nu]*psiexp(-j + nu,1) + vEG2G2G2G2cADBH*b[kappa]*b[nu]*psiexp(kappa - nu,2)*psiexp(-j + nu,1) + vEG2G2G2G2cAWBD*b[kappa]*b[nu]*psiexp(j - k,1)*psiexp(-j + nu,2) + vEG2G2G2G2cAWBDCH*b[kappa]*b[nu]*psiexp(j - k,1)*psiexp(kappa - nu,3)*psiexp(-j + nu,2) + vEG2G2G2G2cAWAHBD*b[kappa]*b[nu]*psiexp(j - k + kappa - nu,1)*psiexp(-j + nu,2) + vEG2G2G2G2cAWAD*b[kappa]*b[nu]*psiexp(-k + nu,1) + vEG2G2G2G2cAWADBH*b[kappa]*b[nu]*psiexp(kappa - nu,2)*psiexp(-k + nu,1);
				}
			}
			MM(j, k) += sumInuIkappa;
		}
	}
}

for (int j = 2; j <= -2 + q; j++)
{
	for (int k = 1; k <= -1 + j; k++)
	{
		{
			scalar sumInuIkappa = 0;
			for (int nu = 2 + j; nu <= q; nu++)
			{
				for (int kappa = 1 + j; kappa <= -1 + nu; kappa++)
				{
					sumInuIkappa += vEG2G2G2G2c*b[kappa]*b[nu] + vEG2G2G2G2cAW*b[kappa]*b[nu]*psiexp(j - k,1) + vEG2G2G2G2cAD*b[kappa]*b[nu]*psiexp(-j + kappa,1) + vEG2G2G2G2cAWBD*b[kappa]*b[nu]*psiexp(j - k,1)*psiexp(-j + kappa,2) + vEG2G2G2G2cAWAD*b[kappa]*b[nu]*psiexp(-k + kappa,1) + vEG2G2G2G2cADAH*b[kappa]*b[nu]*psiexp(-j + nu,1) + vEG2G2G2G2cAWBDBH*b[kappa]*b[nu]*psiexp(j - k,1)*psiexp(-j + nu,2) + vEG2G2G2G2cAWADAH*b[kappa]*b[nu]*psiexp(-k + nu,1) + vEG2G2G2G2cAH*b[kappa]*b[nu]*psiexp(-kappa + nu,1) + vEG2G2G2G2cADBH*b[kappa]*b[nu]*psiexp(-j + kappa,1)*psiexp(-kappa + nu,2) + vEG2G2G2G2cAWADBH*b[kappa]*b[nu]*psiexp(-k + kappa,1)*psiexp(-kappa + nu,2) + vEG2G2G2G2cAWBDCH*b[kappa]*b[nu]*psiexp(j - k,1)*psiexp(-j + kappa,2)*psiexp(-kappa + nu,3) + vEG2G2G2G2cAWAH*b[kappa]*b[nu]*psiexp(j - k - kappa + nu,1) + vEG2G2G2G2cAWAHBD*b[kappa]*b[nu]*psiexp(-j + kappa,2)*psiexp(j - k - kappa + nu,1);
				}
			}
			MM(j, k) += sumInuIkappa;
		}
	}
}

for (int j = 3; j <= -1 + q; j++)
{
	for (int k = 1; k <= -2 + j; k++)
	{
		{
			scalar sumInuIkappa = 0;
			for (int nu = 1 + j; nu <= q; nu++)
			{
				for (int kappa = 1 + k; kappa <= -1 + j; kappa++)
				{
					sumInuIkappa += vEG2G2G2G2c*b[kappa]*b[nu] + vEG2G2G2G2cAWAD*b[kappa]*b[nu]*psiexp(j - k,1) + vEG2G2G2G2cAD*b[kappa]*b[nu]*psiexp(j - kappa,1) + vEG2G2G2G2cAW*b[kappa]*b[nu]*psiexp(-k + kappa,1) + vEG2G2G2G2cAWBD*b[kappa]*b[nu]*psiexp(j - kappa,2)*psiexp(-k + kappa,1) + vEG2G2G2G2cAH*b[kappa]*b[nu]*psiexp(-j + nu,1) + vEG2G2G2G2cAWADBH*b[kappa]*b[nu]*psiexp(j - k,1)*psiexp(-j + nu,2) + vEG2G2G2G2cADBH*b[kappa]*b[nu]*psiexp(j - kappa,1)*psiexp(-j + nu,2) + vEG2G2G2G2cAWBDCH*b[kappa]*b[nu]*psiexp(j - kappa,2)*psiexp(-k + kappa,1)*psiexp(-j + nu,3) + vEG2G2G2G2cAWADAH*b[kappa]*b[nu]*psiexp(-k + nu,1) + vEG2G2G2G2cADAH*b[kappa]*b[nu]*psiexp(-kappa + nu,1) + vEG2G2G2G2cAWBDBH*b[kappa]*b[nu]*psiexp(-k + kappa,1)*psiexp(-kappa + nu,2) + vEG2G2G2G2cAWAH*b[kappa]*b[nu]*psiexp(-j - k + kappa + nu,1) + vEG2G2G2G2cAWAHBD*b[kappa]*b[nu]*psiexp(j - kappa,2)*psiexp(-j - k + kappa + nu,1);
				}
			}
			MM(j, k) += sumInuIkappa;
		}
	}
}

for (int j = 2; j <= -1 + q; j++)
{
	for (int k = 1; k <= -1 + j; k++)
	{
		{
			scalar sumInuIkappa = 0;
			for (int nu = 1 + j; nu <= q; nu++)
			{
				for (int kappa = 0; kappa <= -1 + k; kappa++)
				{
					sumInuIkappa += vEG2G2G2G2c*b[kappa]*b[nu] + vEG2G2G2G2cAD*b[kappa]*b[nu]*psiexp(j - k,1) + vEG2G2G2G2cAWAD*b[kappa]*b[nu]*psiexp(j - kappa,1) + vEG2G2G2G2cAW*b[kappa]*b[nu]*psiexp(k - kappa,1) + vEG2G2G2G2cAWBD*b[kappa]*b[nu]*psiexp(j - k,2)*psiexp(k - kappa,1) + vEG2G2G2G2cAH*b[kappa]*b[nu]*psiexp(-j + nu,1) + vEG2G2G2G2cADBH*b[kappa]*b[nu]*psiexp(j - k,1)*psiexp(-j + nu,2) + vEG2G2G2G2cAWADBH*b[kappa]*b[nu]*psiexp(j - kappa,1)*psiexp(-j + nu,2) + vEG2G2G2G2cAWBDCH*b[kappa]*b[nu]*psiexp(j - k,2)*psiexp(k - kappa,1)*psiexp(-j + nu,3) + vEG2G2G2G2cADAH*b[kappa]*b[nu]*psiexp(-k + nu,1) + vEG2G2G2G2cAWBDBH*b[kappa]*b[nu]*psiexp(k - kappa,1)*psiexp(-k + nu,2) + vEG2G2G2G2cAWADAH*b[kappa]*b[nu]*psiexp(-kappa + nu,1) + vEG2G2G2G2cAWAH*b[kappa]*b[nu]*psiexp(-j + k - kappa + nu,1) + vEG2G2G2G2cAWAHBD*b[kappa]*b[nu]*psiexp(j - k,2)*psiexp(-j + k - kappa + nu,1);
				}
			}
			MM(j, k) += sumInuIkappa;
		}
	}
}

for (int j = 2; j <= -1 + q; j++)
{
	for (int k = 1; k <= -1 + j; k++)
	{
		{
			scalar sumInuIkappa = 0;
			for (int nu = 1 + j; nu <= q; nu++)
			{
				for (int kappa = nu; kappa <= nu; kappa++)
				{
					sumInuIkappa += vEG2G2G4c*Power2(b[nu]) + vEG2G2G4cAW*Power2(b[nu])*psiexp(j - k,1) + vEG2G2G4cAD*Power2(b[nu])*psiexp(-j + nu,1) + vEG2G2G4cAWBD*Power2(b[nu])*psiexp(j - k,1)*psiexp(-j + nu,2) + vEG2G2G4cAWAD*Power2(b[nu])*psiexp(-k + nu,1);
				}
			}
			MM(j, k) += sumInuIkappa;
		}
	}
}

for (int j = 2; j <= -1 + q; j++)
{
	for (int k = 1; k <= -1 + j; k++)
	{
		{
			scalar sumInuIkappa = 0;
			for (int nu = 1 + j; nu <= q; nu++)
			{
				for (int kappa = j; kappa <= j; kappa++)
				{
					sumInuIkappa += vEG2G4G2c*b[kappa]*b[nu] + vEG2G4G2cAW*b[kappa]*b[nu]*psiexp(-k + kappa,1) + vEG2G4G2cAWAH*b[kappa]*b[nu]*psiexp(-k + nu,1) + vEG2G4G2cAH*b[kappa]*b[nu]*psiexp(-kappa + nu,1) + vEG2G4G2cBH*b[kappa]*b[nu]*psiexp(-kappa + nu,2) + vEG2G4G2cAWBH*b[kappa]*b[nu]*psiexp(-k + kappa,1)*psiexp(-kappa + nu,2) + vEG2G4G2cAWCH*b[kappa]*b[nu]*psiexp(-k + kappa,1)*psiexp(-kappa + nu,3);
				}
			}
			MM(j, k) += sumInuIkappa;
		}
	}
}

for (int j = 2; j <= -1 + q; j++)
{
	for (int k = 1; k <= -1 + j; k++)
	{
		{
			scalar sumInuIkappa = 0;
			for (int nu = 1 + j; nu <= q; nu++)
			{
				for (int kappa = k; kappa <= k; kappa++)
				{
					sumInuIkappa += vEG4G2G2c*b[kappa]*b[nu] + vEG4G2G2cAD*b[kappa]*b[nu]*psiexp(j - kappa,1) + vEG4G2G2cBD*b[kappa]*b[nu]*psiexp(j - kappa,2) + vEG4G2G2cAH*b[kappa]*b[nu]*psiexp(-j + nu,1) + vEG4G2G2cAHBD*b[kappa]*b[nu]*psiexp(j - kappa,2)*psiexp(-j + nu,1) + vEG4G2G2cADBH*b[kappa]*b[nu]*psiexp(j - kappa,1)*psiexp(-j + nu,2) + vEG4G2G2cBDCH*b[kappa]*b[nu]*psiexp(j - kappa,2)*psiexp(-j + nu,3) + vEG4G2G2cADAH*b[kappa]*b[nu]*psiexp(-kappa + nu,1) + vEG4G2G2cBDBH*b[kappa]*b[nu]*psiexp(-kappa + nu,2);
				}
			}
			MM(j, k) += sumInuIkappa;
		}
	}
}

for (int j = 3; j <= -1 + q; j++)
{
	for (int k = 1; k <= -2 + j; k++)
	{
		{
			scalar sumInuIkappa = 0;
			for (int nu = 1 + k; nu <= -1 + j; nu++)
			{
				for (int kappa = 1 + j; kappa <= q; kappa++)
				{
					sumInuIkappa += vEG2G2G2G2c*b[kappa]*b[nu] + vEG2G2G2G2cAWAD*b[kappa]*b[nu]*psiexp(j - k,1) + vEG2G2G2G2cAH*b[kappa]*b[nu]*psiexp(-j + kappa,1) + vEG2G2G2G2cAWADBH*b[kappa]*b[nu]*psiexp(j - k,1)*psiexp(-j + kappa,2) + vEG2G2G2G2cAWADAH*b[kappa]*b[nu]*psiexp(-k + kappa,1) + vEG2G2G2G2cAD*b[kappa]*b[nu]*psiexp(j - nu,1) + vEG2G2G2G2cADBH*b[kappa]*b[nu]*psiexp(-j + kappa,2)*psiexp(j - nu,1) + vEG2G2G2G2cADAH*b[kappa]*b[nu]*psiexp(kappa - nu,1) + vEG2G2G2G2cAW*b[kappa]*b[nu]*psiexp(-k + nu,1) + vEG2G2G2G2cAWBD*b[kappa]*b[nu]*psiexp(j - nu,2)*psiexp(-k + nu,1) + vEG2G2G2G2cAWBDCH*b[kappa]*b[nu]*psiexp(-j + kappa,3)*psiexp(j - nu,2)*psiexp(-k + nu,1) + vEG2G2G2G2cAWBDBH*b[kappa]*b[nu]*psiexp(kappa - nu,2)*psiexp(-k + nu,1) + vEG2G2G2G2cAWAH*b[kappa]*b[nu]*psiexp(-j - k + kappa + nu,1) + vEG2G2G2G2cAWAHBD*b[kappa]*b[nu]*psiexp(j - nu,2)*psiexp(-j - k + kappa + nu,1);
				}
			}
			MM(j, k) += sumInuIkappa;
		}
	}
}

for (int j = 4; j <= q; j++)
{
	for (int k = 1; k <= -3 + j; k++)
	{
		{
			scalar sumInuIkappa = 0;
			for (int nu = 1 + k; nu <= -2 + j; nu++)
			{
				for (int kappa = 1 + nu; kappa <= -1 + j; kappa++)
				{
					sumInuIkappa += vEG2G2G2G2c*b[kappa]*b[nu] + vEG2G2G2G2cAWADAH*b[kappa]*b[nu]*psiexp(j - k,1) + vEG2G2G2G2cAH*b[kappa]*b[nu]*psiexp(j - kappa,1) + vEG2G2G2G2cAWAD*b[kappa]*b[nu]*psiexp(-k + kappa,1) + vEG2G2G2G2cAWADBH*b[kappa]*b[nu]*psiexp(j - kappa,2)*psiexp(-k + kappa,1) + vEG2G2G2G2cADAH*b[kappa]*b[nu]*psiexp(j - nu,1) + vEG2G2G2G2cAD*b[kappa]*b[nu]*psiexp(kappa - nu,1) + vEG2G2G2G2cADBH*b[kappa]*b[nu]*psiexp(j - kappa,2)*psiexp(kappa - nu,1) + vEG2G2G2G2cAW*b[kappa]*b[nu]*psiexp(-k + nu,1) + vEG2G2G2G2cAWBDBH*b[kappa]*b[nu]*psiexp(j - nu,2)*psiexp(-k + nu,1) + vEG2G2G2G2cAWBD*b[kappa]*b[nu]*psiexp(kappa - nu,2)*psiexp(-k + nu,1) + vEG2G2G2G2cAWBDCH*b[kappa]*b[nu]*psiexp(j - kappa,3)*psiexp(kappa - nu,2)*psiexp(-k + nu,1) + vEG2G2G2G2cAWAH*b[kappa]*b[nu]*psiexp(j - k - kappa + nu,1) + vEG2G2G2G2cAWAHBD*b[kappa]*b[nu]*psiexp(kappa - nu,2)*psiexp(j - k - kappa + nu,1);
				}
			}
			MM(j, k) += sumInuIkappa;
		}
	}
}

for (int j = 4; j <= q; j++)
{
	for (int k = 1; k <= -3 + j; k++)
	{
		{
			scalar sumInuIkappa = 0;
			for (int nu = 2 + k; nu <= -1 + j; nu++)
			{
				for (int kappa = 1 + k; kappa <= -1 + nu; kappa++)
				{
					sumInuIkappa += vEG2G2G2G2c*b[kappa]*b[nu] + vEG2G2G2G2cAWADAH*b[kappa]*b[nu]*psiexp(j - k,1) + vEG2G2G2G2cADAH*b[kappa]*b[nu]*psiexp(j - kappa,1) + vEG2G2G2G2cAW*b[kappa]*b[nu]*psiexp(-k + kappa,1) + vEG2G2G2G2cAWBDBH*b[kappa]*b[nu]*psiexp(j - kappa,2)*psiexp(-k + kappa,1) + vEG2G2G2G2cAH*b[kappa]*b[nu]*psiexp(j - nu,1) + vEG2G2G2G2cAWAH*b[kappa]*b[nu]*psiexp(j - k + kappa - nu,1) + vEG2G2G2G2cAWAD*b[kappa]*b[nu]*psiexp(-k + nu,1) + vEG2G2G2G2cAWADBH*b[kappa]*b[nu]*psiexp(j - nu,2)*psiexp(-k + nu,1) + vEG2G2G2G2cAD*b[kappa]*b[nu]*psiexp(-kappa + nu,1) + vEG2G2G2G2cADBH*b[kappa]*b[nu]*psiexp(j - nu,2)*psiexp(-kappa + nu,1) + vEG2G2G2G2cAWBD*b[kappa]*b[nu]*psiexp(-k + kappa,1)*psiexp(-kappa + nu,2) + vEG2G2G2G2cAWBDCH*b[kappa]*b[nu]*psiexp(-k + kappa,1)*psiexp(j - nu,3)*psiexp(-kappa + nu,2) + vEG2G2G2G2cAWAHBD*b[kappa]*b[nu]*psiexp(j - k + kappa - nu,1)*psiexp(-kappa + nu,2);
				}
			}
			MM(j, k) += sumInuIkappa;
		}
	}
}

for (int j = 3; j <= q; j++)
{
	for (int k = 1; k <= -2 + j; k++)
	{
		{
			scalar sumInuIkappa = 0;
			for (int nu = 1 + k; nu <= -1 + j; nu++)
			{
				for (int kappa = 0; kappa <= -1 + k; kappa++)
				{
					sumInuIkappa += vEG2G2G2G2c*b[kappa]*b[nu] + vEG2G2G2G2cADAH*b[kappa]*b[nu]*psiexp(j - k,1) + vEG2G2G2G2cAWADAH*b[kappa]*b[nu]*psiexp(j - kappa,1) + vEG2G2G2G2cAW*b[kappa]*b[nu]*psiexp(k - kappa,1) + vEG2G2G2G2cAWBDBH*b[kappa]*b[nu]*psiexp(j - k,2)*psiexp(k - kappa,1) + vEG2G2G2G2cAH*b[kappa]*b[nu]*psiexp(j - nu,1) + vEG2G2G2G2cAWAH*b[kappa]*b[nu]*psiexp(j + k - kappa - nu,1) + vEG2G2G2G2cAD*b[kappa]*b[nu]*psiexp(-k + nu,1) + vEG2G2G2G2cADBH*b[kappa]*b[nu]*psiexp(j - nu,2)*psiexp(-k + nu,1) + vEG2G2G2G2cAWBD*b[kappa]*b[nu]*psiexp(k - kappa,1)*psiexp(-k + nu,2) + vEG2G2G2G2cAWBDCH*b[kappa]*b[nu]*psiexp(k - kappa,1)*psiexp(j - nu,3)*psiexp(-k + nu,2) + vEG2G2G2G2cAWAHBD*b[kappa]*b[nu]*psiexp(j + k - kappa - nu,1)*psiexp(-k + nu,2) + vEG2G2G2G2cAWAD*b[kappa]*b[nu]*psiexp(-kappa + nu,1) + vEG2G2G2G2cAWADBH*b[kappa]*b[nu]*psiexp(j - nu,2)*psiexp(-kappa + nu,1);
				}
			}
			MM(j, k) += sumInuIkappa;
		}
	}
}

for (int j = 3; j <= q; j++)
{
	for (int k = 1; k <= -2 + j; k++)
	{
		{
			scalar sumInuIkappa = 0;
			for (int nu = 1 + k; nu <= -1 + j; nu++)
			{
				for (int kappa = j; kappa <= j; kappa++)
				{
					sumInuIkappa += vEG2G2G4c*b[kappa]*b[nu] + vEG2G2G4cAWAD*b[kappa]*b[nu]*psiexp(-k + kappa,1) + vEG2G2G4cAD*b[kappa]*b[nu]*psiexp(kappa - nu,1) + vEG2G2G4cAW*b[kappa]*b[nu]*psiexp(-k + nu,1) + vEG2G2G4cAWBD*b[kappa]*b[nu]*psiexp(kappa - nu,2)*psiexp(-k + nu,1);
				}
			}
			MM(j, k) += sumInuIkappa;
		}
	}
}

for (int j = 3; j <= q; j++)
{
	for (int k = 1; k <= -2 + j; k++)
	{
		{
			scalar sumInuIkappa = 0;
			for (int nu = 1 + k; nu <= -1 + j; nu++)
			{
				for (int kappa = nu; kappa <= nu; kappa++)
				{
					sumInuIkappa += vEG2G4G2c*Power2(b[nu]) + vEG2G4G2cAWAH*Power2(b[nu])*psiexp(j - k,1) + vEG2G4G2cAH*Power2(b[nu])*psiexp(j - nu,1) + vEG2G4G2cBH*Power2(b[nu])*psiexp(j - nu,2) + vEG2G4G2cAW*Power2(b[nu])*psiexp(-k + nu,1) + vEG2G4G2cAWBH*Power2(b[nu])*psiexp(j - nu,2)*psiexp(-k + nu,1) + vEG2G4G2cAWCH*Power2(b[nu])*psiexp(j - nu,3)*psiexp(-k + nu,1);
				}
			}
			MM(j, k) += sumInuIkappa;
		}
	}
}

for (int j = 3; j <= q; j++)
{
	for (int k = 1; k <= -2 + j; k++)
	{
		{
			scalar sumInuIkappa = 0;
			for (int nu = 1 + k; nu <= -1 + j; nu++)
			{
				for (int kappa = k; kappa <= k; kappa++)
				{
					sumInuIkappa += vEG4G2G2c*b[kappa]*b[nu] + vEG4G2G2cADAH*b[kappa]*b[nu]*psiexp(j - kappa,1) + vEG4G2G2cBDBH*b[kappa]*b[nu]*psiexp(j - kappa,2) + vEG4G2G2cAH*b[kappa]*b[nu]*psiexp(j - nu,1) + vEG4G2G2cAD*b[kappa]*b[nu]*psiexp(-kappa + nu,1) + vEG4G2G2cADBH*b[kappa]*b[nu]*psiexp(j - nu,2)*psiexp(-kappa + nu,1) + vEG4G2G2cBD*b[kappa]*b[nu]*psiexp(-kappa + nu,2) + vEG4G2G2cAHBD*b[kappa]*b[nu]*psiexp(j - nu,1)*psiexp(-kappa + nu,2) + vEG4G2G2cBDCH*b[kappa]*b[nu]*psiexp(j - nu,3)*psiexp(-kappa + nu,2);
				}
			}
			MM(j, k) += sumInuIkappa;
		}
	}
}

for (int j = 2; j <= -1 + q; j++)
{
	for (int k = 1; k <= -1 + j; k++)
	{
		{
			scalar sumInuIkappa = 0;
			for (int nu = 0; nu <= -1 + k; nu++)
			{
				for (int kappa = 1 + j; kappa <= q; kappa++)
				{
					sumInuIkappa += vEG2G2G2G2c*b[kappa]*b[nu] + vEG2G2G2G2cAD*b[kappa]*b[nu]*psiexp(j - k,1) + vEG2G2G2G2cAH*b[kappa]*b[nu]*psiexp(-j + kappa,1) + vEG2G2G2G2cADBH*b[kappa]*b[nu]*psiexp(j - k,1)*psiexp(-j + kappa,2) + vEG2G2G2G2cADAH*b[kappa]*b[nu]*psiexp(-k + kappa,1) + vEG2G2G2G2cAWAD*b[kappa]*b[nu]*psiexp(j - nu,1) + vEG2G2G2G2cAWADBH*b[kappa]*b[nu]*psiexp(-j + kappa,2)*psiexp(j - nu,1) + vEG2G2G2G2cAW*b[kappa]*b[nu]*psiexp(k - nu,1) + vEG2G2G2G2cAWBD*b[kappa]*b[nu]*psiexp(j - k,2)*psiexp(k - nu,1) + vEG2G2G2G2cAWBDCH*b[kappa]*b[nu]*psiexp(j - k,2)*psiexp(-j + kappa,3)*psiexp(k - nu,1) + vEG2G2G2G2cAWBDBH*b[kappa]*b[nu]*psiexp(-k + kappa,2)*psiexp(k - nu,1) + vEG2G2G2G2cAWADAH*b[kappa]*b[nu]*psiexp(kappa - nu,1) + vEG2G2G2G2cAWAH*b[kappa]*b[nu]*psiexp(-j + k + kappa - nu,1) + vEG2G2G2G2cAWAHBD*b[kappa]*b[nu]*psiexp(j - k,2)*psiexp(-j + k + kappa - nu,1);
				}
			}
			MM(j, k) += sumInuIkappa;
		}
	}
}

for (int j = 3; j <= q; j++)
{
	for (int k = 1; k <= -2 + j; k++)
	{
		{
			scalar sumInuIkappa = 0;
			for (int nu = 0; nu <= -1 + k; nu++)
			{
				for (int kappa = 1 + k; kappa <= -1 + j; kappa++)
				{
					sumInuIkappa += vEG2G2G2G2c*b[kappa]*b[nu] + vEG2G2G2G2cADAH*b[kappa]*b[nu]*psiexp(j - k,1) + vEG2G2G2G2cAH*b[kappa]*b[nu]*psiexp(j - kappa,1) + vEG2G2G2G2cAD*b[kappa]*b[nu]*psiexp(-k + kappa,1) + vEG2G2G2G2cADBH*b[kappa]*b[nu]*psiexp(j - kappa,2)*psiexp(-k + kappa,1) + vEG2G2G2G2cAWADAH*b[kappa]*b[nu]*psiexp(j - nu,1) + vEG2G2G2G2cAW*b[kappa]*b[nu]*psiexp(k - nu,1) + vEG2G2G2G2cAWBDBH*b[kappa]*b[nu]*psiexp(j - k,2)*psiexp(k - nu,1) + vEG2G2G2G2cAWBD*b[kappa]*b[nu]*psiexp(-k + kappa,2)*psiexp(k - nu,1) + vEG2G2G2G2cAWBDCH*b[kappa]*b[nu]*psiexp(j - kappa,3)*psiexp(-k + kappa,2)*psiexp(k - nu,1) + vEG2G2G2G2cAWAH*b[kappa]*b[nu]*psiexp(j + k - kappa - nu,1) + vEG2G2G2G2cAWAHBD*b[kappa]*b[nu]*psiexp(-k + kappa,2)*psiexp(j + k - kappa - nu,1) + vEG2G2G2G2cAWAD*b[kappa]*b[nu]*psiexp(kappa - nu,1) + vEG2G2G2G2cAWADBH*b[kappa]*b[nu]*psiexp(j - kappa,2)*psiexp(kappa - nu,1);
				}
			}
			MM(j, k) += sumInuIkappa;
		}
	}
}

for (int j = 3; j <= q; j++)
{
	for (int k = 2; k <= -1 + j; k++)
	{
		{
			scalar sumInuIkappa = 0;
			for (int nu = 0; nu <= -2 + k; nu++)
			{
				for (int kappa = 1 + nu; kappa <= -1 + k; kappa++)
				{
					sumInuIkappa += vEG2G2G2G2c*b[kappa]*b[nu] + vEG2G2G2G2cAH*b[kappa]*b[nu]*psiexp(j - k,1) + vEG2G2G2G2cADAH*b[kappa]*b[nu]*psiexp(j - kappa,1) + vEG2G2G2G2cAD*b[kappa]*b[nu]*psiexp(k - kappa,1) + vEG2G2G2G2cADBH*b[kappa]*b[nu]*psiexp(j - k,2)*psiexp(k - kappa,1) + vEG2G2G2G2cAWADAH*b[kappa]*b[nu]*psiexp(j - nu,1) + vEG2G2G2G2cAWAD*b[kappa]*b[nu]*psiexp(k - nu,1) + vEG2G2G2G2cAWADBH*b[kappa]*b[nu]*psiexp(j - k,2)*psiexp(k - nu,1) + vEG2G2G2G2cAW*b[kappa]*b[nu]*psiexp(kappa - nu,1) + vEG2G2G2G2cAWBDBH*b[kappa]*b[nu]*psiexp(j - kappa,2)*psiexp(kappa - nu,1) + vEG2G2G2G2cAWBD*b[kappa]*b[nu]*psiexp(k - kappa,2)*psiexp(kappa - nu,1) + vEG2G2G2G2cAWBDCH*b[kappa]*b[nu]*psiexp(j - k,3)*psiexp(k - kappa,2)*psiexp(kappa - nu,1) + vEG2G2G2G2cAWAH*b[kappa]*b[nu]*psiexp(j - k + kappa - nu,1) + vEG2G2G2G2cAWAHBD*b[kappa]*b[nu]*psiexp(k - kappa,2)*psiexp(j - k + kappa - nu,1);
				}
			}
			MM(j, k) += sumInuIkappa;
		}
	}
}

for (int j = 3; j <= q; j++)
{
	for (int k = 2; k <= -1 + j; k++)
	{
		{
			scalar sumInuIkappa = 0;
			for (int nu = 1; nu <= -1 + k; nu++)
			{
				for (int kappa = 0; kappa <= -1 + nu; kappa++)
				{
					sumInuIkappa += vEG2G2G2G2c*b[kappa]*b[nu] + vEG2G2G2G2cAH*b[kappa]*b[nu]*psiexp(j - k,1) + vEG2G2G2G2cAWADAH*b[kappa]*b[nu]*psiexp(j - kappa,1) + vEG2G2G2G2cAWAD*b[kappa]*b[nu]*psiexp(k - kappa,1) + vEG2G2G2G2cAWADBH*b[kappa]*b[nu]*psiexp(j - k,2)*psiexp(k - kappa,1) + vEG2G2G2G2cADAH*b[kappa]*b[nu]*psiexp(j - nu,1) + vEG2G2G2G2cAD*b[kappa]*b[nu]*psiexp(k - nu,1) + vEG2G2G2G2cADBH*b[kappa]*b[nu]*psiexp(j - k,2)*psiexp(k - nu,1) + vEG2G2G2G2cAW*b[kappa]*b[nu]*psiexp(-kappa + nu,1) + vEG2G2G2G2cAWBDBH*b[kappa]*b[nu]*psiexp(j - nu,2)*psiexp(-kappa + nu,1) + vEG2G2G2G2cAWBD*b[kappa]*b[nu]*psiexp(k - nu,2)*psiexp(-kappa + nu,1) + vEG2G2G2G2cAWBDCH*b[kappa]*b[nu]*psiexp(j - k,3)*psiexp(k - nu,2)*psiexp(-kappa + nu,1) + vEG2G2G2G2cAWAH*b[kappa]*b[nu]*psiexp(j - k - kappa + nu,1) + vEG2G2G2G2cAWAHBD*b[kappa]*b[nu]*psiexp(k - nu,2)*psiexp(j - k - kappa + nu,1);
				}
			}
			MM(j, k) += sumInuIkappa;
		}
	}
}

for (int j = 2; j <= q; j++)
{
	for (int k = 1; k <= -1 + j; k++)
	{
		{
			scalar sumInuIkappa = 0;
			for (int nu = 0; nu <= -1 + k; nu++)
			{
				for (int kappa = j; kappa <= j; kappa++)
				{
					sumInuIkappa += vEG2G2G4c*b[kappa]*b[nu] + vEG2G2G4cAD*b[kappa]*b[nu]*psiexp(-k + kappa,1) + vEG2G2G4cAW*b[kappa]*b[nu]*psiexp(k - nu,1) + vEG2G2G4cAWBD*b[kappa]*b[nu]*psiexp(-k + kappa,2)*psiexp(k - nu,1) + vEG2G2G4cAWAD*b[kappa]*b[nu]*psiexp(kappa - nu,1);
				}
			}
			MM(j, k) += sumInuIkappa;
		}
	}
}

for (int j = 2; j <= q; j++)
{
	for (int k = 1; k <= -1 + j; k++)
	{
		{
			scalar sumInuIkappa = 0;
			for (int nu = 0; nu <= -1 + k; nu++)
			{
				for (int kappa = k; kappa <= k; kappa++)
				{
					sumInuIkappa += vEG2G4G2c*b[kappa]*b[nu] + vEG2G4G2cAH*b[kappa]*b[nu]*psiexp(j - kappa,1) + vEG2G4G2cBH*b[kappa]*b[nu]*psiexp(j - kappa,2) + vEG2G4G2cAWAH*b[kappa]*b[nu]*psiexp(j - nu,1) + vEG2G4G2cAW*b[kappa]*b[nu]*psiexp(kappa - nu,1) + vEG2G4G2cAWBH*b[kappa]*b[nu]*psiexp(j - kappa,2)*psiexp(kappa - nu,1) + vEG2G4G2cAWCH*b[kappa]*b[nu]*psiexp(j - kappa,3)*psiexp(kappa - nu,1);
				}
			}
			MM(j, k) += sumInuIkappa;
		}
	}
}

for (int j = 2; j <= q; j++)
{
	for (int k = 1; k <= -1 + j; k++)
	{
		{
			scalar sumInuIkappa = 0;
			for (int nu = 0; nu <= -1 + k; nu++)
			{
				for (int kappa = nu; kappa <= nu; kappa++)
				{
					sumInuIkappa += vEG4G2G2c*Power2(b[nu]) + vEG4G2G2cAH*Power2(b[nu])*psiexp(j - k,1) + vEG4G2G2cADAH*Power2(b[nu])*psiexp(j - nu,1) + vEG4G2G2cBDBH*Power2(b[nu])*psiexp(j - nu,2) + vEG4G2G2cAD*Power2(b[nu])*psiexp(k - nu,1) + vEG4G2G2cADBH*Power2(b[nu])*psiexp(j - k,2)*psiexp(k - nu,1) + vEG4G2G2cBD*Power2(b[nu])*psiexp(k - nu,2) + vEG4G2G2cAHBD*Power2(b[nu])*psiexp(j - k,1)*psiexp(k - nu,2) + vEG4G2G2cBDCH*Power2(b[nu])*psiexp(j - k,3)*psiexp(k - nu,2);
				}
			}
			MM(j, k) += sumInuIkappa;
		}
	}
}

for (int j = 2; j <= -1 + q; j++)
{
	for (int k = 1; k <= -1 + j; k++)
	{
		{
			scalar sumInuIkappa = 0;
			for (int nu = j; nu <= j; nu++)
			{
				for (int kappa = 1 + j; kappa <= q; kappa++)
				{
					sumInuIkappa += vEG2G4G2c*b[kappa]*b[nu] + vEG2G4G2cAWAH*b[kappa]*b[nu]*psiexp(-k + kappa,1) + vEG2G4G2cAH*b[kappa]*b[nu]*psiexp(kappa - nu,1) + vEG2G4G2cBH*b[kappa]*b[nu]*psiexp(kappa - nu,2) + vEG2G4G2cAW*b[kappa]*b[nu]*psiexp(-k + nu,1) + vEG2G4G2cAWBH*b[kappa]*b[nu]*psiexp(kappa - nu,2)*psiexp(-k + nu,1) + vEG2G4G2cAWCH*b[kappa]*b[nu]*psiexp(kappa - nu,3)*psiexp(-k + nu,1);
				}
			}
			MM(j, k) += sumInuIkappa;
		}
	}
}

for (int j = 3; j <= q; j++)
{
	for (int k = 1; k <= -2 + j; k++)
	{
		{
			scalar sumInuIkappa = 0;
			for (int nu = j; nu <= j; nu++)
			{
				for (int kappa = 1 + k; kappa <= -1 + j; kappa++)
				{
					sumInuIkappa += vEG2G2G4c*b[kappa]*b[nu] + vEG2G2G4cAW*b[kappa]*b[nu]*psiexp(-k + kappa,1) + vEG2G2G4cAWAD*b[kappa]*b[nu]*psiexp(-k + nu,1) + vEG2G2G4cAD*b[kappa]*b[nu]*psiexp(-kappa + nu,1) + vEG2G2G4cAWBD*b[kappa]*b[nu]*psiexp(-k + kappa,1)*psiexp(-kappa + nu,2);
				}
			}
			MM(j, k) += sumInuIkappa;
		}
	}
}

for (int j = 2; j <= q; j++)
{
	for (int k = 1; k <= -1 + j; k++)
	{
		{
			scalar sumInuIkappa = 0;
			for (int nu = j; nu <= j; nu++)
			{
				for (int kappa = 0; kappa <= -1 + k; kappa++)
				{
					sumInuIkappa += vEG2G2G4c*b[kappa]*b[nu] + vEG2G2G4cAW*b[kappa]*b[nu]*psiexp(k - kappa,1) + vEG2G2G4cAD*b[kappa]*b[nu]*psiexp(-k + nu,1) + vEG2G2G4cAWBD*b[kappa]*b[nu]*psiexp(k - kappa,1)*psiexp(-k + nu,2) + vEG2G2G4cAWAD*b[kappa]*b[nu]*psiexp(-kappa + nu,1);
				}
			}
			MM(j, k) += sumInuIkappa;
		}
	}
}

for (int j = 2; j <= q; j++)
{
	for (int k = 1; k <= -1 + j; k++)
	{
		{
			scalar sumInuIkappa = 0;
			for (int nu = j; nu <= j; nu++)
			{
				for (int kappa = j; kappa <= j; kappa++)
				{
					sumInuIkappa += vEG2G6c*Power2(b[nu]) + vEG2G6cAW*Power2(b[nu])*psiexp(-k + nu,1);
				}
			}
			MM(j, k) += sumInuIkappa;
		}
	}
}

for (int j = 2; j <= q; j++)
{
	for (int k = 1; k <= -1 + j; k++)
	{
		{
			scalar sumInuIkappa = 0;
			for (int nu = j; nu <= j; nu++)
			{
				for (int kappa = k; kappa <= k; kappa++)
				{
					sumInuIkappa += vEG4G4c*b[kappa]*b[nu] + vEG4G4cAD*b[kappa]*b[nu]*psiexp(-kappa + nu,1) + vEG4G4cBD*b[kappa]*b[nu]*psiexp(-kappa + nu,2);
				}
			}
			MM(j, k) += sumInuIkappa;
		}
	}
}

for (int j = 2; j <= -1 + q; j++)
{
	for (int k = 1; k <= -1 + j; k++)
	{
		{
			scalar sumInuIkappa = 0;
			for (int nu = k; nu <= k; nu++)
			{
				for (int kappa = 1 + j; kappa <= q; kappa++)
				{
					sumInuIkappa += vEG4G2G2c*b[kappa]*b[nu] + vEG4G2G2cAH*b[kappa]*b[nu]*psiexp(-j + kappa,1) + vEG4G2G2cAD*b[kappa]*b[nu]*psiexp(j - nu,1) + vEG4G2G2cADBH*b[kappa]*b[nu]*psiexp(-j + kappa,2)*psiexp(j - nu,1) + vEG4G2G2cBD*b[kappa]*b[nu]*psiexp(j - nu,2) + vEG4G2G2cAHBD*b[kappa]*b[nu]*psiexp(-j + kappa,1)*psiexp(j - nu,2) + vEG4G2G2cBDCH*b[kappa]*b[nu]*psiexp(-j + kappa,3)*psiexp(j - nu,2) + vEG4G2G2cADAH*b[kappa]*b[nu]*psiexp(kappa - nu,1) + vEG4G2G2cBDBH*b[kappa]*b[nu]*psiexp(kappa - nu,2);
				}
			}
			MM(j, k) += sumInuIkappa;
		}
	}
}

for (int j = 3; j <= q; j++)
{
	for (int k = 1; k <= -2 + j; k++)
	{
		{
			scalar sumInuIkappa = 0;
			for (int nu = k; nu <= k; nu++)
			{
				for (int kappa = 1 + k; kappa <= -1 + j; kappa++)
				{
					sumInuIkappa += vEG4G2G2c*b[kappa]*b[nu] + vEG4G2G2cAH*b[kappa]*b[nu]*psiexp(j - kappa,1) + vEG4G2G2cADAH*b[kappa]*b[nu]*psiexp(j - nu,1) + vEG4G2G2cBDBH*b[kappa]*b[nu]*psiexp(j - nu,2) + vEG4G2G2cAD*b[kappa]*b[nu]*psiexp(kappa - nu,1) + vEG4G2G2cADBH*b[kappa]*b[nu]*psiexp(j - kappa,2)*psiexp(kappa - nu,1) + vEG4G2G2cBD*b[kappa]*b[nu]*psiexp(kappa - nu,2) + vEG4G2G2cAHBD*b[kappa]*b[nu]*psiexp(j - kappa,1)*psiexp(kappa - nu,2) + vEG4G2G2cBDCH*b[kappa]*b[nu]*psiexp(j - kappa,3)*psiexp(kappa - nu,2);
				}
			}
			MM(j, k) += sumInuIkappa;
		}
	}
}

for (int j = 2; j <= q; j++)
{
	for (int k = 1; k <= -1 + j; k++)
	{
		{
			scalar sumInuIkappa = 0;
			for (int nu = k; nu <= k; nu++)
			{
				for (int kappa = 0; kappa <= -1 + k; kappa++)
				{
					sumInuIkappa += vEG2G4G2c*b[kappa]*b[nu] + vEG2G4G2cAWAH*b[kappa]*b[nu]*psiexp(j - kappa,1) + vEG2G4G2cAH*b[kappa]*b[nu]*psiexp(j - nu,1) + vEG2G4G2cBH*b[kappa]*b[nu]*psiexp(j - nu,2) + vEG2G4G2cAW*b[kappa]*b[nu]*psiexp(-kappa + nu,1) + vEG2G4G2cAWBH*b[kappa]*b[nu]*psiexp(j - nu,2)*psiexp(-kappa + nu,1) + vEG2G4G2cAWCH*b[kappa]*b[nu]*psiexp(j - nu,3)*psiexp(-kappa + nu,1);
				}
			}
			MM(j, k) += sumInuIkappa;
		}
	}
}

for (int j = 2; j <= q; j++)
{
	for (int k = 1; k <= -1 + j; k++)
	{
		{
			scalar sumInuIkappa = 0;
			for (int nu = k; nu <= k; nu++)
			{
				for (int kappa = j; kappa <= j; kappa++)
				{
					sumInuIkappa += vEG4G4c*b[kappa]*b[nu] + vEG4G4cAD*b[kappa]*b[nu]*psiexp(kappa - nu,1) + vEG4G4cBD*b[kappa]*b[nu]*psiexp(kappa - nu,2);
				}
			}
			MM(j, k) += sumInuIkappa;
		}
	}
}

for (int j = 2; j <= q; j++)
{
	for (int k = 1; k <= -1 + j; k++)
	{
		{
			scalar sumInuIkappa = 0;
			for (int nu = k; nu <= k; nu++)
			{
				for (int kappa = k; kappa <= k; kappa++)
				{
					sumInuIkappa += vEG6G2c*Power2(b[nu]) + vEG6G2cAH*Power2(b[nu])*psiexp(j - nu,1) + vEG6G2cBH*Power2(b[nu])*psiexp(j - nu,2) + vEG6G2cCH*Power2(b[nu])*psiexp(j - nu,3);
				}
			}
			MM(j, k) += sumInuIkappa;
		}
	}
}

for (int j = 1; j <= -2 + q; j++)
{
	for (int k = j; k <= j; k++)
	{
		{
			scalar sumInuIkappa = 0;
			for (int nu = 1 + j; nu <= -1 + q; nu++)
			{
				for (int kappa = 1 + nu; kappa <= q; kappa++)
				{
					sumInuIkappa += vEG4G2G2c*b[kappa]*b[nu] + vEG4G2G2cADAH*b[kappa]*b[nu]*psiexp(-k + kappa,1) + vEG4G2G2cBDBH*b[kappa]*b[nu]*psiexp(-k + kappa,2) + vEG4G2G2cAH*b[kappa]*b[nu]*psiexp(kappa - nu,1) + vEG4G2G2cAD*b[kappa]*b[nu]*psiexp(-k + nu,1) + vEG4G2G2cADBH*b[kappa]*b[nu]*psiexp(kappa - nu,2)*psiexp(-k + nu,1) + vEG4G2G2cBD*b[kappa]*b[nu]*psiexp(-k + nu,2) + vEG4G2G2cAHBD*b[kappa]*b[nu]*psiexp(kappa - nu,1)*psiexp(-k + nu,2) + vEG4G2G2cBDCH*b[kappa]*b[nu]*psiexp(kappa - nu,3)*psiexp(-k + nu,2);
				}
			}
			MM(j, k) += sumInuIkappa;
		}
	}
}

for (int j = 1; j <= -2 + q; j++)
{
	for (int k = j; k <= j; k++)
	{
		{
			scalar sumInuIkappa = 0;
			for (int nu = 2 + j; nu <= q; nu++)
			{
				for (int kappa = 1 + j; kappa <= -1 + nu; kappa++)
				{
					sumInuIkappa += vEG4G2G2c*b[kappa]*b[nu] + vEG4G2G2cAD*b[kappa]*b[nu]*psiexp(-k + kappa,1) + vEG4G2G2cBD*b[kappa]*b[nu]*psiexp(-k + kappa,2) + vEG4G2G2cADAH*b[kappa]*b[nu]*psiexp(-k + nu,1) + vEG4G2G2cBDBH*b[kappa]*b[nu]*psiexp(-k + nu,2) + vEG4G2G2cAH*b[kappa]*b[nu]*psiexp(-kappa + nu,1) + vEG4G2G2cAHBD*b[kappa]*b[nu]*psiexp(-k + kappa,2)*psiexp(-kappa + nu,1) + vEG4G2G2cADBH*b[kappa]*b[nu]*psiexp(-k + kappa,1)*psiexp(-kappa + nu,2) + vEG4G2G2cBDCH*b[kappa]*b[nu]*psiexp(-k + kappa,2)*psiexp(-kappa + nu,3);
				}
			}
			MM(j, k) += sumInuIkappa;
		}
	}
}

for (int j = 1; j <= -1 + q; j++)
{
	for (int k = j; k <= j; k++)
	{
		{
			scalar sumInuIkappa = 0;
			for (int nu = 1 + j; nu <= q; nu++)
			{
				for (int kappa = 0; kappa <= -1 + j; kappa++)
				{
					sumInuIkappa += vEG2G4G2c*b[kappa]*b[nu] + vEG2G4G2cAW*b[kappa]*b[nu]*psiexp(k - kappa,1) + vEG2G4G2cAH*b[kappa]*b[nu]*psiexp(-k + nu,1) + vEG2G4G2cBH*b[kappa]*b[nu]*psiexp(-k + nu,2) + vEG2G4G2cAWBH*b[kappa]*b[nu]*psiexp(k - kappa,1)*psiexp(-k + nu,2) + vEG2G4G2cAWCH*b[kappa]*b[nu]*psiexp(k - kappa,1)*psiexp(-k + nu,3) + vEG2G4G2cAWAH*b[kappa]*b[nu]*psiexp(-kappa + nu,1);
				}
			}
			MM(j, k) += sumInuIkappa;
		}
	}
}

for (int j = 1; j <= -1 + q; j++)
{
	for (int k = j; k <= j; k++)
	{
		{
			scalar sumInuIkappa = 0;
			for (int nu = 1 + j; nu <= q; nu++)
			{
				for (int kappa = nu; kappa <= nu; kappa++)
				{
					sumInuIkappa += vEG4G4c*Power2(b[nu]) + vEG4G4cAD*Power2(b[nu])*psiexp(-k + nu,1) + vEG4G4cBD*Power2(b[nu])*psiexp(-k + nu,2);
				}
			}
			MM(j, k) += sumInuIkappa;
		}
	}
}

for (int j = 1; j <= -1 + q; j++)
{
	for (int k = j; k <= j; k++)
	{
		{
			scalar sumInuIkappa = 0;
			for (int nu = 1 + j; nu <= q; nu++)
			{
				for (int kappa = j; kappa <= j; kappa++)
				{
					sumInuIkappa += vEG6G2c*b[kappa]*b[nu] + vEG6G2cAH*b[kappa]*b[nu]*psiexp(-kappa + nu,1) + vEG6G2cBH*b[kappa]*b[nu]*psiexp(-kappa + nu,2) + vEG6G2cCH*b[kappa]*b[nu]*psiexp(-kappa + nu,3);
				}
			}
			MM(j, k) += sumInuIkappa;
		}
	}
}

for (int j = 1; j <= -1 + q; j++)
{
	for (int k = j; k <= j; k++)
	{
		{
			scalar sumInuIkappa = 0;
			for (int nu = 0; nu <= -1 + j; nu++)
			{
				for (int kappa = 1 + j; kappa <= q; kappa++)
				{
					sumInuIkappa += vEG2G4G2c*b[kappa]*b[nu] + vEG2G4G2cAH*b[kappa]*b[nu]*psiexp(-k + kappa,1) + vEG2G4G2cBH*b[kappa]*b[nu]*psiexp(-k + kappa,2) + vEG2G4G2cAW*b[kappa]*b[nu]*psiexp(k - nu,1) + vEG2G4G2cAWBH*b[kappa]*b[nu]*psiexp(-k + kappa,2)*psiexp(k - nu,1) + vEG2G4G2cAWCH*b[kappa]*b[nu]*psiexp(-k + kappa,3)*psiexp(k - nu,1) + vEG2G4G2cAWAH*b[kappa]*b[nu]*psiexp(kappa - nu,1);
				}
			}
			MM(j, k) += sumInuIkappa;
		}
	}
}

for (int j = 2; j <= q; j++)
{
	for (int k = j; k <= j; k++)
	{
		{
			scalar sumInuIkappa = 0;
			for (int nu = 0; nu <= -2 + j; nu++)
			{
				for (int kappa = 1 + nu; kappa <= -1 + j; kappa++)
				{
					sumInuIkappa += vEG2G2G4c*b[kappa]*b[nu] + vEG2G2G4cAD*b[kappa]*b[nu]*psiexp(k - kappa,1) + vEG2G2G4cAWAD*b[kappa]*b[nu]*psiexp(k - nu,1) + vEG2G2G4cAW*b[kappa]*b[nu]*psiexp(kappa - nu,1) + vEG2G2G4cAWBD*b[kappa]*b[nu]*psiexp(k - kappa,2)*psiexp(kappa - nu,1);
				}
			}
			MM(j, k) += sumInuIkappa;
		}
	}
}

for (int j = 2; j <= q; j++)
{
	for (int k = j; k <= j; k++)
	{
		{
			scalar sumInuIkappa = 0;
			for (int nu = 1; nu <= -1 + j; nu++)
			{
				for (int kappa = 0; kappa <= -1 + nu; kappa++)
				{
					sumInuIkappa += vEG2G2G4c*b[kappa]*b[nu] + vEG2G2G4cAWAD*b[kappa]*b[nu]*psiexp(k - kappa,1) + vEG2G2G4cAD*b[kappa]*b[nu]*psiexp(k - nu,1) + vEG2G2G4cAW*b[kappa]*b[nu]*psiexp(-kappa + nu,1) + vEG2G2G4cAWBD*b[kappa]*b[nu]*psiexp(k - nu,2)*psiexp(-kappa + nu,1);
				}
			}
			MM(j, k) += sumInuIkappa;
		}
	}
}

for (int j = 1; j <= q; j++)
{
	for (int k = j; k <= j; k++)
	{
		{
			scalar sumInuIkappa = 0;
			for (int nu = 0; nu <= -1 + j; nu++)
			{
				for (int kappa = j; kappa <= j; kappa++)
				{
					sumInuIkappa += vEG2G6c*b[kappa]*b[nu] + vEG2G6cAW*b[kappa]*b[nu]*psiexp(kappa - nu,1);
				}
			}
			MM(j, k) += sumInuIkappa;
		}
	}
}

for (int j = 1; j <= q; j++)
{
	for (int k = j; k <= j; k++)
	{
		{
			scalar sumInuIkappa = 0;
			for (int nu = 0; nu <= -1 + j; nu++)
			{
				for (int kappa = nu; kappa <= nu; kappa++)
				{
					sumInuIkappa += vEG4G4c*Power2(b[nu]) + vEG4G4cAD*Power2(b[nu])*psiexp(k - nu,1) + vEG4G4cBD*Power2(b[nu])*psiexp(k - nu,2);
				}
			}
			MM(j, k) += sumInuIkappa;
		}
	}
}

for (int j = 1; j <= -1 + q; j++)
{
	for (int k = j; k <= j; k++)
	{
		{
			scalar sumInuIkappa = 0;
			for (int nu = j; nu <= j; nu++)
			{
				for (int kappa = 1 + j; kappa <= q; kappa++)
				{
					sumInuIkappa += vEG6G2c*b[kappa]*b[nu] + vEG6G2cAH*b[kappa]*b[nu]*psiexp(kappa - nu,1) + vEG6G2cBH*b[kappa]*b[nu]*psiexp(kappa - nu,2) + vEG6G2cCH*b[kappa]*b[nu]*psiexp(kappa - nu,3);
				}
			}
			MM(j, k) += sumInuIkappa;
		}
	}
}

for (int j = 1; j <= q; j++)
{
	for (int k = j; k <= j; k++)
	{
		{
			scalar sumInuIkappa = 0;
			for (int nu = j; nu <= j; nu++)
			{
				for (int kappa = 0; kappa <= -1 + j; kappa++)
				{
					sumInuIkappa += vEG2G6c*b[kappa]*b[nu] + vEG2G6cAW*b[kappa]*b[nu]*psiexp(-kappa + nu,1);
				}
			}
			MM(j, k) += sumInuIkappa;
		}
	}
}

for (int j = 1; j <= q; j++)
{
	for (int k = j; k <= j; k++)
	{
		{
			scalar sumInuIkappa = 0;
			for (int nu = j; nu <= j; nu++)
			{
				for (int kappa = j; kappa <= j; kappa++)
				{
					sumInuIkappa += vEG8c*Power2(b[nu]);
				}
			}
			MM(j, k) += sumInuIkappa;
		}
	}
}

free(psiexp_array);
#undef Psi
#undef psiexp
