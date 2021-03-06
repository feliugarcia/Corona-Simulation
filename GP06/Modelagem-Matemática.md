Modelagem matemática do projeto
===============================

Caracterização do algoritmo de pandemia
---------------------------------------

Deseja-se criar um modelo espacial (ou de superfície ) que represente o
processo de infecção, em uma população de tamanho n, por um determinado
um determinado agente viral, no caso, coronavirus.

Para construir um processo que represente a propagação da infecção em
uma população, num intervalo de tempo, devemos considerar
características da contaminação, da população em risco e das políticas
adotadas para gestão da crise. Os elementos destacados em cada tópico
abaixo podem ajudar no delineamento da modelagem.

1 . Como características da contaminação, levantamos os seguintes
pontos:

    a . Como as pessoas são infectadas;
    b . Com que velocidade a infecção se espalha na população
    c . Qual é o período de encubação (tempo até a apresentação dos sinais);
    d . Qual é o período a partir do qual um indivíduo passa a gerar contaminação (ambiente ou outros entes;
    e . Qual o impacto dos indivíduos assintomáticos no processo de contaminação;

2 . Como características da população em riscos consideramos:

    a . Mecanismos/canais de contaminação mais frequentes;
    b . Fatores de risco na população, comorbidades;
    c . Capacidade dos sistemas de saúde para atender e tratar os enfermos;

3 . Capacidade do Estado em gerenciar a crise pandêmica

    a . Disponibilidade de instalações, equipamentos e profissionais de para atender os enfermos;
    b . Gabinete integrado de gestão de crise que possibilite acionamento dos diversos equipamentos de governo para prover soluções para tratamento da crise;
    c . Disponibilidade de informações integradas e periódicas sobre o tema que possibilitem a melhor tomada de decisão em qualquer instante da crise; Exemplo de informação necessária é a mensuração da capacidade do sistema para atender à demanda de enfermos em cada instante no tempo.

Certamente podem ser levantados outros pontos potencialmente relevantes
para a modelagem que eventualmente não foram considerados nesta breve
caracterização do problema, no entanto, representam ponto de partida
construção de um protótipo pode trazer um valor para o produto entregue.

Delimitação do problema
-----------------------

Avaliar o impacto de medidas administrativas no controle da propagação
da infecção na população. São consideradas as seguintes medidas:

    a ) isolamento social; 
    b ) minimização das interações sociais; e
    c ) combinação das estratégias “a” e “b”.

### Localidade objeto do experimento:

Para realização do experimento deste trabalho foi selecionada a cidade
de Santana do Jacaré, localizada em Minas Gerais. Com população de 4.613
habitantes, que efetua seus deslocamentos em trechos curtos, e com
limites que determinam um polígono que lembra um retângulo, essa cidade
oferece características adequadas para o esforço de processamento
computacional inicialmente disponível para o projeto.

Descrição do experimento:
-------------------------

Em razão da indisponibilidade de dados sobre deslocamentos de pessoas em
comunidades, há necessidade de propor um modelo de deslocamento que
possa minimamente representar o movimento de transeuntes em uma pequena
cidade, revelando suas interações, em termos de exposição ao vírus, com
o espaço e com indivíduos dos quais se aproximou. Tendo a localização do
indivíduo em cada instante no tempo, aplica-se os procedimentos de
indução e propagação da infecção e recuperação da infecção. Entende-se
que os passos abaixo descritos podem oferecer insights do algoritmo de
pandemia.

1 . Fixar locais de residência dos indivíduos;

2 . Constituir núcleos familiares de mesmo endereço para os indivíduos
da população;

3 . Estabelecer a localização do inicial da população no instante t0,
{(x, y ) pertencente à região da cidade}

4 . Estabelece mecanismo aleatório para determinar a localização da
população no instante t+1.

5 . Selecionar, entre os indivíduos da população, o portador inicial do
vírus (paciente zero).

6 . Estabelecer a superfície de contaminação que seguirá a trilha do
paciente zero.

7 . Induzir a contaminação em indivíduos que estejam nas proximidades do
traço do paciente zero.

8 . Fazer com que os contaminados adotem, depois de algum tempo, os
mesmos comportamentos de contaminação do paciente zero.

9 . Estabelecer mecanismos de cura/recuperação

10 . Estabelecer cenários de redução de interação e de isolamento social

11 . Mensurar o nível de contaminação da população em cada instante no
tempo, considerando a aplicação isolada ou combinada dos cenários de

    a ) redução da iteração;e

    b ) isolamento social.
