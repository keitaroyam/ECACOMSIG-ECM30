library(ggplot2)

# https://kohske.wordpress.com/2010/12/29/using-jet-colormap-in-r-and-ggplot2/
jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))

d<-read.table("zeta.dat",h=T)
p <- ggplot(d, aes(x=x,y=y,fill=zeta)) +geom_tile()
p <- p +scale_fill_gradientn(colours = jet.colors(10))
p <- p +coord_fixed() +theme_void() +theme(legend.position=c(1.1, .8))

ggsave("zeta_map.png", p, w=9, h=7)
