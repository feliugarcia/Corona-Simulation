#############################################################
#   Visualização Lista 2 - Explorando Visualizacoes.Rmd     #
#############################################################

    Este arquivo visar trazer uma explanação do código criado
pelo grupo GP0402--Animacoes.
    As alterações foram realizadas no Chunk:
         ### Animação - Atuação do Grupo GP0402.
         
    Este Chunk se inicia na linha 157, e pode alterar conforme
a evolução do código.

#############################################################
#           Atividades a Realizadas até 03/04/2020          #
#############################################################
Feito:
1) Criação de Animação e arquivo "./Animacao/Deslocamento_Populacao.gif"
2) Método de renderização utilizado: gifski_renderer()
3) Inclusão do pacote "gifski" no arquivo Funcoes_R.R, no p_load do pacman
4) Definição da resolução da imagem renderizada de 1024x600
5) Inclusão de título para a renderização: "Animação da movimentação da população"
6) Definição do tamanho dos pontos para size=.5
7) Diferenciação se o indivíduo está ou não em casa
    - Triângulo: indivíduo em casa. dados.plot$casa = 1
    - Círculo: indivíduo se movimentando. dados.plot$casa = 0
8) Identificação indivíduo infectado
    - Escala de 0 à 1, 0 = indivíduo não infectado; 1 = Infectado
    - Estala a partir de: dados.plot$status
    - Escala baseado em cor. color=status
9) Sorteio de 199 indivíduos da população
10) Frames por segundo = fps = 10
11) Inclusão de Mapa da cidade de Santana do Jacaré/MG
  

A fazer:
1) Inclusão de contador de tempo (frame)
  Possíveis métodos (à verificar e podendo haver outros):
    - lable, grobTree, textGrob, geom_text, annotate
Ref: http://www.sthda.com/english/wiki/ggplot2-texts-add-text-annotations-to-a-graph-in-r-software
