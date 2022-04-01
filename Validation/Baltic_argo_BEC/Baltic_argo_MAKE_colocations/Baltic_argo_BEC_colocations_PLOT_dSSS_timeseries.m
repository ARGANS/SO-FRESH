




PRODS = {'SSS_{NS}','SSS_{NM}','SSS_{gl}'};
 
figure
h1 = plot(dSSS_gl_irange_mean,'k'); hold on
h2 = plot(dSSS_nm_irange_mean,'r','linewidth',2);
h3 = plot(dSSS_ns_irange_mean,'b','linewidth',2);

lg = legend([h1(1),h2(1),h3(1)],PRODS,'fontsize',18,'location','SouthEast');


ylim([-2 4])