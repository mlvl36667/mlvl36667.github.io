document.addEventListener('DOMContentLoaded', () => {
  const imgs = Array.from(document.querySelectorAll('.gallery img'));
  const lb = document.getElementById('lightbox');
  const lbImg = lb.querySelector('img');
  const counter = lb.querySelector('.lightbox-counter');
  let idx = 0;

  function show(i) {
    idx = (i + imgs.length) % imgs.length;
    lbImg.src = imgs[idx].src;
    counter.textContent = (idx + 1) + ' / ' + imgs.length;
  }

  document.querySelector('.gallery').addEventListener('click', e => {
    if (e.target.tagName === 'IMG') { idx = imgs.indexOf(e.target); show(idx); lb.classList.add('active'); }
  });

  lb.querySelector('.lightbox-close').addEventListener('click', () => lb.classList.remove('active'));
  lb.querySelector('.lightbox-prev').addEventListener('click', e => { e.stopPropagation(); show(idx - 1); });
  lb.querySelector('.lightbox-next').addEventListener('click', e => { e.stopPropagation(); show(idx + 1); });
  lb.addEventListener('click', e => { if (e.target === lb) lb.classList.remove('active'); });

  document.addEventListener('keydown', e => {
    if (!lb.classList.contains('active')) return;
    if (e.key === 'ArrowLeft') show(idx - 1);
    else if (e.key === 'ArrowRight') show(idx + 1);
    else if (e.key === 'Escape') lb.classList.remove('active');
  });
});
