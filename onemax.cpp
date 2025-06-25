#include <iostream>
#include <random>
#include <vector>

using std::cout;
using std::endl;
using std::random_device;
using std::mt19937;
using std::uniform_int_distribution;
using std::uniform_real_distribution;
using std::vector;

random_device rd;
mt19937 engine(rd());

struct individual
{
	vector<int> genes;
	double fitness;
	double weight;

	individual(double weight)
	{
		this -> weight = weight;
	}

	bool operator > (const individual &other)
	{
		return this -> fitness * this -> weight > other.fitness * other.weight;
	}

	bool operator < (const individual &other)
	{
		return this -> fitness * this -> weight < other.fitness * other.weight;
	}
	
};

vector<individual> tournament(vector<individual> &pop, int tournsize)
{
        vector<individual> aspirants;
        int pop_size = pop.size();
        uniform_int_distribution dist(0, pop_size - 1);

        for(int i = 0; i < pop_size; i++){
                int one = dist(engine), another_one = 0;
                for(int g = 0; g < tournsize - 1; g++){
                        do{
                                another_one = dist(engine);
                        }while(one == another_one);
                        one = pop[one] > pop[another_one] ? one : another_one;
                }
                aspirants.push_back(pop[one]);
        }
        return aspirants;
}

void crossover(vector<individual> &aspirants, double CXPB)
{
        int pop_size = aspirants.size(), len_gen = aspirants[0].genes.size();
        uniform_real_distribution pb(0.0, 1.0);
        uniform_int_distribution cross_dot(1, len_gen - 2);

        for(int i = 0; i < pop_size - 1; i += 2){
                if(pb(engine) < CXPB){
                        int dot = cross_dot(engine);
                        for(int g = 0; g < dot; g++){
                                std::swap(aspirants[i].genes[g], aspirants[i + 1].genes[g]);
                        }
                }
        }
}

void mutation(vector<individual> &offspring, double MTPB, double MGPB)
{
        int pop_size = offspring.size(), len_gen = offspring[0].genes.size();
        uniform_real_distribution pb(0.0, 1.0);
        uniform_int_distribution gene(0, len_gen - 1);

        for(int i = 0; i < pop_size; i++){
                if(pb(engine) < MTPB){
                        for(int g = 0; g < len_gen; g++){
                                if(pb(engine) < MGPB){
                                        offspring[i].genes[g] = !offspring[i].genes[g];
                                }
                        }
                }
        }
}

double stat(vector<individual> &pop)
{
        int pop_size = pop.size();
        double min = pop[0].fitness, max = pop[0].fitness, avg = pop[0].fitness;

        for(int i = 1; i < pop_size; i++){
                min = pop[i].fitness < min ? pop[i].fitness : min;
                max = pop[i].fitness > max ? pop[i].fitness : max;
                avg += pop[i].fitness;
        }
        avg /= pop_size;
        cout << min << '\t' << max << '\t' << avg << endl;
	return max;
}

double evalute(individual &individ)
{
	double sum = 0;

	for(int g = 0; g < individ.genes.size(); g++){
		sum += individ.genes[g];
	}

	return sum;
}

void compute_fitness(vector<individual> &population)
{
	int pop_size = population.size();

	for(int i = 0; i < pop_size; i++){
		population[i].fitness= evalute(population[i]);
	}
}

int main(void)
{
        int eras = 500;
        int len_genes = 100, pop_size = 2000, tournsize = 3;
        double CXPB = 0.7, MTPB = 0.2, MGPB = 0.001;

        uniform_int_distribution zero_or_one(0, 1);

        vector<individual> population;

        for(int i = 0; i < pop_size; i++){
                individual individ(1.0);
                for(int g = 0; g < len_genes; g++){
                        individ.genes.push_back(zero_or_one(engine));
                }
                individ.fitness = evalute(individ);
                population.push_back(individ);
        }

        cout << "era\tmin\tmax\tavg" << endl;
        cout << 0 << '\t';
        stat(population);

	for(int g = 1; g < eras + 1; g++){
		population = tournament(population, tournsize);
		crossover(population, CXPB);
		mutation(population, MTPB, MGPB);
		compute_fitness(population);
		cout << g << '\t';
		double max = stat(population);
		if(max == len_genes)
			break;
	}

        return 0;
}
