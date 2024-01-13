%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% This code is for the manuscript                                       %%%%%%
%%%%%% 'Effects of environmental feedback on species with finite population' %%%%%%
%%%%%% submitted to iScience.                                                %%%%%%
%%%%%% Code author: Jia-Xu Han                                               %%%%%%
%%%%%% 2024.1.12                                                             %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc
species_size=2; %Number of species
enviroment_size=species_size;%Number of environments (resource)
population_size=50; %Population size for the species of concern
parameter_idx=50; %Number of parts to divide the parameter into
individual_number=100; %Total individual number for all species
parameter_c=ones(species_size,1);    %Parameter c in equation (4) for all species
parameter_d=ones(species_size,1);    %Parameter d in equation (4) for all species
parameter_selection=0.9; %Selection intensity
parameter_a_collection=linspace(0,1,parameter_idx)'; %The collection of parameter a in equation (5) for environment
parameter_c(1)=0.5; %The collection of parameter c in equation (4) for the species of concern
result_sum=zeros(parameter_idx,parameter_idx); %Matrix for recording results
repeat=10000; %Number of repetitions per simulation
T_max=1001;  %Maximum time step to run in each simulation
dt=0.01; %Delta t in equation (5)
result_sum_env=zeros(parameter_idx,T_max-1);
for idx3=1:parameter_idx
    parameter_a=ones(species_size,1)*parameter_a_collection(idx3); %Parameter a in equation (5) for environment
    result_record_env=zeros(repeat,T_max-1); %Vector for temporary recording of environmental state
    for idx2=1:repeat
        population_size_simulation=[population_size;(individual_number-population_size)]; %Initialize the population size
        environment_simulation=ones(enviroment_size,1)/(enviroment_size); %Initialize the environment
        T=1;  %Initialize simulation time-record
        while(population_size_simulation(1)~=0 & population_size_simulation(1)~=individual_number & T<T_max)
            density_population=population_size_simulation/individual_number; %Calculate the population density of each species
            change_environment=parameter_a.*environment_simulation.*(1-environment_simulation).*(density_population-flip(density_population)); %Calculate the environment change in each time step
            payoff_species=parameter_c.*environment_simulation+parameter_d; %Calculate the pay-off of each species
            fitness_species=1-parameter_selection+parameter_selection*payoff_species;  %Calculate the fitness of each species
            probability_birth=fitness_species.*population_size_simulation/sum(fitness_species.*population_size_simulation);  %Calculate the birth rate of each species
            choose_birth=sum(sign(sign(rand-cumsum(probability_birth))+1))+1;  %Randomly select one species to produce an offspring
            choose_death=sum(sign(sign(rand-cumsum(density_population))+1))+1;  %Randomly select one individual to die
            population_size_simulation(choose_birth)=population_size_simulation(choose_birth)+1; %Updata population state
            population_size_simulation(choose_death)=population_size_simulation(choose_death)-1; %Updata population state
            environment_simulation=environment_simulation+change_environment*dt; %Updata environment state
            result_record_env(idx2,T)=environment_simulation(1); 
            T=T+1;
        end
    end
    result_sum_env(idx3,:)=mean(sign(sign(1.5*result_record_env-1)+1));
end
figure(1)
imagesc([1:T_max],(parameter_a_collection),result_sum_env)
xlabel('Time')
ylabel('a')
colorbar
