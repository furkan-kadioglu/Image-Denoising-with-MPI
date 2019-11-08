/*
Student Name: Furkan Kadıoğlu
Student Number: 2015400051
Compile Status: Compiling
Program Status: Working
*/


#include <iostream>
#include <vector>
#include <mpi.h>
#include <fstream>
#include <math.h>
#include <random>


using namespace std;


/*
    GLOBAL VARIABLE FOR GENERALITY
*/
//shows dimension of mesh 
int dimension = 2;

//size of image is fixed
int row_size = 200;
int col_size = 200;

//column & row size of processors
int per_col;
int per_row;

//number of processors per column || row
int pro_col;
int pro_row;

//program values
double beta,pi,_gamma;
int num_of_step = 500000;



/*
    FUNCTION PART
    'random_generator' creates random number at given interval 
    'get_who' returns rank of processor which has given index
*/


//generates random number at the interval which is wanted 
double random_generator(double low, double high)
{   
    random_device rd;  //Will be used to obtain a seed for the random number engine
    mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    uniform_real_distribution<> dis(low, high);

    return dis(gen);
}



//returns which processor has given index 
int get_who(int row, int col)
{  
    return (row / per_row) * pro_col + (col / per_col) + 1;
}



/*  
    MAIN FUNCTION
    global variables are calculated @main
    what master and slaves will do are defined @main 
*/
    


int main(int argc,char* const argv[])
{
     // Initialize the MPI environment
    MPI_Init(NULL, NULL);
    
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);



    //initiliaze beta & pi values
    beta = atof(argv[3]);
    pi = atof(argv[4]);    
    _gamma = 0.5 * log((1 - pi) / pi);

    //number of processor per column || row 
    pro_col = pow(sqrt(world_size-1),(dimension-1));
    pro_row = pow((world_size - 1),1.0/dimension);

    //keeps how many pixels are kept @one processor
    per_col = col_size / pro_col;
    per_row = row_size / pro_row;

    
    
    /*
        MASTER
    */
    if(world_rank == 0)
    {
        //reads pixel values and pushes image vector
        ifstream inFile(argv[1]);
        
        /*
            MASTER PART I
            Noised pixel values are distributed to slaves
        */
        // ith processor 
        for(int i=0; i<pro_row; i++)
        {   
            int pix;
            //jth row in ith processor
            for(int j=0; j<per_row; j++)
            {   
                //kth column processor
                for(int k=0; k<pro_col; k++)
                {
                    //temporary memory 
                    vector<int> img;
                    //lth column in kth column procesor
                    for(int l=0; l<per_col; l++)
                    {
                        inFile >> pix;
                        img.push_back(pix);
                    }
                    MPI_Send(&img[0],per_col,MPI_INT,(i*pro_col)+(k+1),0,MPI_COMM_WORLD);
                }
                
                

            }
        }
        inFile.close();

        

        /*
            MASTER PART II
            Denoised pixel values are collected from slaves and 
            are written to output file
            
        */
       
        ofstream outFile(argv[2]);
        int im=0;
        for(int i=0; i<pro_row; i++)
        {
            for(int j=0; j<per_row; j++)
            {
                for(int k=0; k<pro_col; k++)
                {
                    vector<int> out(per_col);
                    
                    MPI_Recv(&out[0],per_col,MPI_INT,(i*pro_col)+(k+1),0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                    
                    for(int l=0; l<per_col; l++)
                    {
                        if((k != pro_col - 1) || (l != per_col - 1))
                            outFile << out[l] << " ";
                        else 
                            outFile << out[l];
                    }
                    
                }
                outFile<<endl;
                
            }
        }
        outFile.close();
        
    }

    /*
        SLAVE PROCESSOR
    */
    else
    {   
        
        /*
            SLAVE PART I
            Noised pixel values are received from master
        */
        

        //keeps noised pixel values 
        vector<vector<int>> img_pro(per_row);
        //gets pixels
        for(int i=0; i<per_row; i++)
        {   
            //temporary memory 
            vector<int> img_col(per_col);
            MPI_Recv(&img_col[0],per_col,MPI_INT,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            img_pro[i] = img_col;
            
        }

        //copied pixel values before iterations
        vector<vector<int>> copy = img_pro;



        /*
            SLAVE PART II
            Finds and keeps neighbors before iterations 
        */
        
        vector<int> neighbors;
        map<int,vector<pair<int,int>>> requested_pixels;
        for(int i=0; i<dimension; i++)
        {   
            
            for(int j=(-1)*i; j<=i; j++)
            {   
             
                if(i != 1 || ((world_rank % pro_col != 0 || j != i) && (world_rank % pro_col != 1 || j != -i)))
                {
                        if(i != 0 || world_rank % pro_col != 0 || dimension == 1)
                        {
                            int n = j + pow(pro_col,(i)) + world_rank;
                            if(n > 0 && n < world_size)
                            {
                                neighbors.push_back(n);
                                requested_pixels[n];
                            }
                                
                        }
                        
                            
                        if(i != 0 || world_rank % pro_col != 1 || dimension == 1)
                        {
                            int n = j - pow(pro_col,(i)) + world_rank;
                            if(n > 0 && n < world_size)
                            {
                                neighbors.push_back(n);
                                requested_pixels[n];
                            }
                                
                        }

                }
                
            }
        }
        sort(neighbors.begin(),neighbors.end());
        
        
        
        /*
            SLAVE PART III
            Iterations are made to denoise
        */
        
       
        
        for(int z=0; z<num_of_step / (world_size - 1); z++)
        {
            
            
            
            /*
                ITERATION PART I
                checks is there any needs from higher priority process 
            */
            
            int next = 0;
            //temporary memory to recevive ciao
            if(neighbors[0] < world_rank)
            {
                
                vector<pair<int,int>> inquire(3);
                for(int i = next; neighbors[i] < world_rank && i < neighbors.size(); i++)
                {  

                    MPI_Recv(&inquire[0],6,MPI_INT,neighbors[i],0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                    
                    if(inquire[0].first != -1)
                    {
                        vector<int> pixel_req;
                        for(auto j: inquire)
                        {
                            if(j.first != -1)
                                pixel_req.push_back(img_pro[j.first][j.second]);
                            else
                                break;
                        }

                        //completes pixel_req size to 3
                        while(pixel_req.size() < 3)
                            pixel_req.push_back(0);


                        MPI_Send(&pixel_req[0],3,MPI_INT,neighbors[i],0,MPI_COMM_WORLD);
                    }
                    next++;
                }
            }
            
            
        


            /*
                ITERATIONS PART II
                checks random pixel is noisy or not  
            */

            //selects pixel to process 
            int row = (int)random_generator(0,per_row);
            int col = (int)random_generator(0,per_col);
            
            //transform to general index
            row = row + (per_row * ((world_rank - 1) / pro_col));
            col = col + (per_col * ((world_rank - 1) % pro_col));
            int pixel = img_pro[row%per_row][col%per_col];
            
            
            
            //collects pixel values to calculate delta_E
            int env = 0;
            
            for(int i = row - 1; i <= row + 1; i++)
            {
                for(int j = col - 1; j <= col + 1; j++)
                {   
                    int i_temp = min(max(i,0),row_size-1);
                    int j_temp = min(max(j,0),col_size-1);
                    //temporary variable
                    int who = get_who(i_temp,j_temp);
                    
                    //transform to local index
                    int r = i_temp % per_row;
                    int c = j_temp % per_col;
                    if(r < 0)
                        r += per_row;
                    if(c < 0) 
                        c += per_col;

                    
                    if(who == world_rank)
                        env += img_pro[r][c];  

                    else 
                        requested_pixels[who].push_back(make_pair(r,c));
                    
                        
                        
                        
                        
                    
                }
            }
            
         
  
            //completes vector size to three
            for(auto &m: requested_pixels)
            {   
                while(m.second.size() < 3)
                    m.second.push_back(make_pair(-1,-1));
            }
            
            
            //requests from neighbors corresponding pixel
            for(auto &i: requested_pixels)
            {   
                MPI_Send(&i.second[0],6,MPI_INT,i.first,0,MPI_COMM_WORLD);
                if(i.second[0].first != -1)
                {
                    vector<int> pixel_values(3);
                    MPI_Recv(&pixel_values[0],3,MPI_INT,i.first,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                    for(auto j:pixel_values)
                        env += j;
                        
                    
                }
                
            }
            for(auto &i : requested_pixels)
                i.second.clear();
            
            env -= pixel;
            
            
            //calculates delta_E
            int multi = copy[row%per_row][col%per_col] * pixel;
            double delta_E = -2*_gamma*multi -2*beta*pixel*env;
            
            //decides flip or not 
            if(log(random_generator(0,1)) < delta_E)
                img_pro[row%per_row][col%per_col] *= -1;
            
            
                
            
            
            
            /*
                ITERATION PART III 
                checks is there any needs from lower priority process 
            */


            
            if(neighbors[neighbors.size() - 1] > world_rank)
            {
          
                vector<pair<int,int>> inquire(3);
                for(int i = next; i < neighbors.size(); i++)
                {
                    MPI_Recv(&inquire[0],6,MPI_INT,neighbors[i],0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                    

                    if(inquire[0].first != -1)
                    {
                        vector<int> pixel_req;
                        for(auto j:inquire)
                        {
                            if(j.first != -1)
                                pixel_req.push_back(img_pro[j.first][j.second]);
                            else
                                break;
                        }

                        //completes pixel_req size to 3
                        while(pixel_req.size() < 3)
                            pixel_req.push_back(0);


                        MPI_Send(&pixel_req[0],3,MPI_INT,neighbors[i],0,MPI_COMM_WORLD);
                    }   
                }
            }
            
            
            
            
        }
        

    
        
        /*
            SLAVE PART IV
            Denoised pixel values are sent to master
        */
        
        for(int i=0; i<per_row; i++)
            MPI_Send(&img_pro[i][0],per_col,MPI_INT,0,0,MPI_COMM_WORLD);
        

    }


    MPI_Finalize();
    
    
    return 0;

}