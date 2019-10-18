fitness=0;
                ispen=false; % si es 1 tiene penalidad
                
                % Chequeo de cruce con el área de cobertura principal:
                % añade inhibidores y resta si estan fuera del área de cobertura out areas:
                for n=1:L % por cada inhibidor
                    ind1=ind(n);
                    aaa=aa(ind1);
                    bbb=bb(ind1);
                    if rot(n)
                        tmp=aaa;
                        aaa=bbb;
                        bbb=tmp;
                    end
                    A0=AA(ind1); % área a cubrir
                    x1=max([x-aaa/2  0]);
                    y1=max([y-bbb/2  0]);
                    x2=min([x+aaa/2  a]);
                    y2=min([y+bbb/2  b]);
                    % x1 - x2,  y1 - y2 son las coordenada del inhibidor
                    if (x1>=x2)||(y1>=y2)
                        A=0; % inhibidor que está adentro del área a inhibir
                    else
                        A=((x2-x1)*(y2-y1)); % inhibidor dentro del área de cobertura
                       
                    end
                    if A<A0 % si no está totalmente adentro del área de cobertura
                        fitness=fitness + A-nac*(A0-A);
                     
                        ispen=true; % penalidad
                    else
                        fitness=fitness + A;
                     
                    end
                end
                
                % por cada par de inhibidores:
                for n1=1:L-1
                    ind1=ind(n1);
                    aaa1=aa(ind1);
                    bbb1=bb(ind1);
                    if rot(n1)
                        tmp=aaa1;
                        aaa1=bbb1;
                        bbb1=tmp;
                    end
                    A1=AA(ind1);
                    x1=x(n1);
                    y1=y(n1); % posición del primer inhibidor
                    for n2=n1+1:L
                        ind2=ind(n2);
                        aaa2=aa(ind2);
                        bbb2=bb(ind2);
                        if rot(n2)
                            tmp=aaa2;
                            aaa2=bbb2;
                            bbb2=tmp;
                        end
                        A2=AA(ind2);
                        x2=x(n2);
                        y2=y(n2); % posición del segundo inhibidor
                        dx=abs(x1-x2);
                        dy=abs(y1-y2); % distancias
                        a12=(aaa1/2+aaa2/2);
                        b12=(bbb1/2+bbb2/2);
                        if (dx<a12)&&(dy<b12) % si cruza
                            ispen=true;
                            Ac=((a12-dx)*(b12-dy)); % area de cruce
                            fitness=fitness-Ac-Ac; % por que el área de n1 y n2 fue añadida completamente
                            fitness=fitness-2*nac*Ac;
                        end

                    end
                end
                
                if ispen
                    fitness=fitness-penalty;
                end