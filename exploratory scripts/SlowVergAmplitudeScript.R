t %>%
  filter(monkey %in% c('Kopachuck','Patos','Pilchuck','QT')) %>%
  group_by(neuron) %>%
  do(preparetoBOOTslow(.)) ->
  t


t %>%
  group_by(neuron) %>%
  summarize(n_slow=length(unique(slowverg)))->
  slow_verg_n_strabismus

ggplot(slow_verg_n_strabismus)+
  geom_histogram(aes(n_slow),binwidth = 10)


t %>%
  group_by(neuron,slowverg) %>%
  summarize(verg_amp=first(verg.amp)) %>%
  rename(move_ID=slowverg) %>%
  separate(neuron,c('monkey','cellnum'),remove=TRUE)->
  slow_verg_amp_strabismus


ggplot(slow_verg_amp_strabismus)+
  geom_histogram(aes(verg_amp),color='white',binwidth = 1)+
  theme_minimal()+
  xlim(-15,15)+
  xlab('Amplitude of Slow Vergence Movement (deg)')+
  labs(caption='Minimum amplitude of accepted movement was 0.5 deg\nOnly strabismic animals included')

write.csv(slow_verg_n_strabismus,'NumberOfSlowVergenceMovements.csv')
write.csv(slow_verg_amp_strabismus,'AmplitudeofSlowVergenceMovements.csv')
