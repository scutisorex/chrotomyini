# Plots for figures 3 and 4
chr.1.1.plot|chr.2.1.plot|chr.3.1.plot|chr.4.1.plot

two.1 <- chr.2.1.plot+theme(axis.title.y=element_blank(),
                            axis.text.y=element_blank(),
                            axis.ticks.y=element_blank())

three.1 <- chr.3.1.plot+theme(axis.title.y=element_blank(),
                            axis.text.y=element_blank(),
                            axis.ticks.y=element_blank())

four.1 <- chr.4.1.plot+theme(axis.title.y=element_blank(),
                              axis.text.y=element_blank(),
                              axis.ticks.y=element_blank())
chr.1.1.plot|two.1|three.1|four.1



two <- chr.2.plot+theme(axis.title.y=element_blank(),
                            axis.text.y=element_blank(),
                            axis.ticks.y=element_blank())

three <- chr.3.plot+theme(axis.title.y=element_blank(),
                              axis.text.y=element_blank(),
                              axis.ticks.y=element_blank())

four <- chr.4.plot+theme(axis.title.y=element_blank(),
                             axis.text.y=element_blank(),
                             axis.ticks.y=element_blank())
chr.1.plot|two|three|four
