using Random

# -- Função Genetic Algorithm --
# Argumentos:
#    p - precisão
#    N - numero de cromossomos
#    zmin - limite inferior
#    zmax - limite superior
#    fun - função a maximizar

function GA(p,N,zmin,zmax,fun,pc=.9,pm=.02,It_max=1000)
    m = 1 + log(2,(zmax-zmin)*10^p)
    m = Int(ceil(m))
    Pop = bitrand(m,N)
    MAX = []
    MED = []
    MIN = []
    it = 0
    Pop_R = zeros(N)

    while true
        for i = 1:N
            Pop_R[i] = bin2real(Pop[:,i],zmin,zmax)
        end
        f = abs.(fun.(Pop_R))
        MAX = [MAX;maximum(fun.(Pop_R))]
        MIN = [MIN;minimum(fun.(Pop_R))]
        MED = [MED;mean(fun.(Pop_R))]
        best_pos = findmax(fun.(Pop_R))[2]
        best = Pop[:,best_pos]
        F = sum(f)
        p = f/F
        q = zeros(N); q[1] = p[1]
        for i = 2:N
            q[i] = q[i-1]+p[i]
        end
        Pop = selecao(Pop,q)
        Pop = crossover(Pop,pc)
        Pop = mutation(Pop,pm)
        Pop[:,1] = best
        if it == It_max
            break
        end
        dens = zeros(100)
        for i = 1:N
            for j = 1:100
                xinf = zmin + j*(zmax-zmin)/100
                xsup = zmin + (j+1)*(zmax-zmin)/100
                if Pop_R[i]>xinf && Pop_R[i]<=xsup
                    dens[j]+=1
                end
            end
        end
        display(maximum(dens))
        if maximum(dens)>=.5*N
            break
        end
        it += 1
    end
    display(Pop)
    p=plot(MAX,m=2,label="Máximo",xlabel="Geração",ylabel="Melhor resultado",title="Gráfico de tendência do AG")
    p=plot!(MIN,m=2,label="Mínimo")
    p=plot!(MED,m=2,label="Média")
    display(p)
    p=plot(fun,xlim=(zmin,zmax),label="Função",title="Função de avaliação de teste")
    p=scatter!(Pop_R,fun.(Pop_R),m=2,label="Geração $it")
    p=scatter!([Pop_R[1]],[fun.(Pop_R)[1]],m=3,label="Máximo encontrado")
    display(p)
end

fun1(x)=21.5 + exp(x)*sin(4pi*x)
fun2(x)=1+exp(-x^2)*cos(36x)
fun3(x)=x*sin(10pi*x)+1

function bin2real(bin,zmin,zmax)
    m = length(bin)
    real = 0
    for i = 1:m
        real = real + bin[i]*2^(i-1)
    end
    real = zmin + real*(zmax-zmin)/(2^m-1)
    return real
end

function selecao(Pop,q)
    V=zeros(size(Pop))
    N = size(Pop,2)
    for i = 1:N
        r = rand()
        if r<=q[1]
            V[:,i] = Pop[:,1]
        else
            for j = 1:(N-1)
                if r>q[j] && r<q[j+1]
                    V[:,i]=Pop[:,j]
                end
            end
        end
    end
    return Bool.(V)
end

function crossover(Pop,pc)
    N = size(Pop,2)
    m = size(Pop,1)
    V= zeros(size(Pop))
    for i =1:(N-1)
        k = rand(1:m)
        r = rand()
        if r <= pc
            V[1:k,i]=Pop[1:k,i]
            V[k+1:end,i]=Pop[k+1:end,i+1]
        end
    end
    V[:,end]=Pop[:,end]
    return Bool.(V)
end

function mutation(Pop,pm)
    N = size(Pop,2)
    m = size(Pop,1)
    for i =1:N
        for j = 1:m
            r = rand()
            if r<=pm
                Pop[j,i] = Pop[j,i]==false
            end
        end
    end
    return Pop
end
