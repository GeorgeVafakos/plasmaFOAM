#include<stdio.h>
#include<math.h>

#define g 9.807

// Declaire functions
void clean_stdin(void);
void clear_screen(void);

// Main program
int main ()
{
    int print_Pr=0, print_Gr=0;
    double Re, Ha, Rm, N, Pr, MagPr, Gr;
    double u, L, nu, rho, B, sigma, mu, alpha, beta, DT;
    char yn;
    
    clear_screen();

    printf("\t\t\t\t--- DIMENSIONLESS NUMBERS CALCULATOR ---\n");
    printf("-----------------------------------------------------------------------------------------------------------------\n");
    printf("This is a program that calculates the dimensionless numbers for the thermal, MHD flow.\n");
    printf("Calculated dimensionless numbers:   Re, Ha, Rm, N\n");
    printf("Optional dimensionless numbers:     Pr, Mag_Pr, Gr\n");
    printf("\n");
    printf("Characteristic velocity:          "); 
    scanf("%lf", &u);
    printf("Characteristic length:            ");
    scanf("%lf", &L);
    printf("Viscosity:                        ");
    scanf("%lf", &nu);
    printf("Density:                          ");
    scanf("%lf", &rho);
    printf("Characteristic Magnetic field:    ");
    scanf("%lf", &B);
    printf("Magnetic permeability:            ");
    scanf("%lf", &mu);
    printf("Electrical conductivity:          ");
    scanf("%lf", &sigma);
    printf("\n");

    clean_stdin();

    while(1)
    {     
        printf("Calculate Prandtl Number? [y/N]   ");
        yn=fgetc(stdin);
        if (yn=='y' || yn=='Y')
        {
            printf("Thermal diffusivity:              ");
            scanf("%lf", &alpha);
            print_Pr = 1;
            clean_stdin();
            break;
        }
        else if (yn=='n' || yn=='N' || yn==0x0A)
            break;
        else
        {
            printf("For YES type 'y' and for NO type 'n' or press enter.\n\n");
            clean_stdin();
        }
    }
    
        
    while(1)
    {
        printf("Calculate Grashof Number? [y/N]   ");
        yn=fgetc(stdin);
        if (yn=='y' || yn=='Y')
        {
            printf("Thermal expansion coefficient:    ");
            scanf("%lf", &beta);
            printf("Temperature difference:           ");
            scanf("%lf", &DT);
            print_Gr = 1;
            break;
        }
        else if (yn=='n' || yn=='N' || yn==0x0A)
            break;
        else
        {
            printf("For YES type 'y' and for NO type 'n' or press enter.\n\n");
            clean_stdin();
        }
    }
    
    // Dimensionless Numbers
    Re = u*L/nu;
    Ha = B*L*sqrt(sigma/(rho*nu));
    Rm = mu*sigma*u*L;
    N  = Ha*Ha/Re;
    
    if (print_Pr == 1)
    {
        Pr = nu/alpha;
        MagPr = nu*mu*sigma;
    }   
    if (print_Gr == 1)
        Gr = g*beta*DT*L*L*L/(nu*nu);
    
    
    printf("\n\n");
    
    printf("--- DIMENSIONLESS NUMBERS ---\n\n");
    
    // Print dimensionless numbers
    printf("The dimensionless numbers for MHD flow are:\n\n");
    printf("Reynolds Number                  Re = %.4f\n",Re);
    printf("Hartmann Number                  Ha = %.1f\n",Ha);
    printf("Magnetic Reynolds Number         Rm = %.4f\n",Rm);
    printf("Interaction Parameter (Stuart)   N  = %.4f\n",N);
    if (print_Pr)
    {
        printf("Prandtl Number                   Pr = %.4f\n",Pr);
        printf("Magnetic Prandtl Number      Pr_mag = %.4f\n",MagPr);
    }
    if (print_Gr)
        printf("Grashof Number                   Gr = %.4f\n",Gr);
    printf("\n\n");
    
}


// Auxiliary functions
void clean_stdin(void)
{
    int c;
    do {
        c = getchar();
    } while (c != '\n' && c != EOF);
}

void clear_screen(void)
{
    printf("\033[H\033[J");
}
